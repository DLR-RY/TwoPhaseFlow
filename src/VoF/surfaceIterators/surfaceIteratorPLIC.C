/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019-2020 DLR
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "surfaceIteratorPLIC.H"
#include "scalarMatrices.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::surfaceIteratorPLIC::surfaceIteratorPLIC
(
    const fvMesh& mesh,
    const scalar tol
)
:
    mesh_(mesh),
    cutCell_(mesh_),
    surfCellTol_(tol)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::label Foam::surfaceIteratorPLIC::vofCutCell
(
    const label celli,
    const scalar alpha1,
    const scalar tol,
    const label maxIter,
    vector normal
)
{
    if (mag(normal) == 0)
    {
        WarningInFunction
            << "normal length is zero in cell: " << celli << nl
            << "try increasing nCorrectors" << endl;

        return sign(alpha1-0.5);
    }

    normal.normalise();

    // Finding cell vertex extrema values
    const labelList& pLabels = mesh_.cellPoints(celli);
    scalarField fvert(pLabels.size());
    forAll(pLabels, pi)
    {
        fvert[pi] = (mesh_.points()[pLabels[pi]] - mesh_.C()[celli]) & (normal);
    }

    labelList order(fvert.size());
    sortedOrder(fvert, order);
    scalar f1 = fvert[order.first()];
    scalar f2 = fvert[order.last()];


    //Handling special case where method is handed an almost full or empty cell
    if (alpha1 < tol)
    {
        return -1;
    }
    else if (1 - alpha1 < tol)
    {
        return 1;
    }

    // Finding the two vertices inbetween which the isovalue giving alpha1 lies
    label L1 = 0;
    label L2 = fvert.size() - 1;
    scalar a1 = 1;
    scalar a2 = 0;
    scalar L3, f3, a3;

    while (L2 - L1 > 1)
    {
        L3 = round(0.5*(L1 + L2));
        f3 = fvert[order[L3]];
        cutCell_.calcSubCell(celli, f3, normal);
        a3 = cutCell_.VolumeOfFluid();
        if (a3 > alpha1)
        {
            L1 = L3; f1 = f3; a1 = a3;
        }
        else if (a3 < alpha1)
        {
            L2 = L3; f2 = f3; a2 = a3;
        }
        else
        {
            return cutCell_.calcSubCell(celli, f3, normal);
        }
    }

    if (mag(f1 - f2) < 10*SMALL)
    {
        return cutCell_.calcSubCell(celli, f1, normal);
    }

    if (mag(a1 - a2) < tol)
    {
        return cutCell_.calcSubCell(celli, 0.5*(f1 + f2), normal);
    }
    // Now we know that a(f) = alpha1 is to be found on the f interval
    // [f1, f2], i.e. alpha1 will be in the interval [a2,a1]


    // Finding coefficients in 3 deg polynomial alpha(f) from 4 solutions

    // Finding 2 additional points on 3 deg polynomial
    f3 = f1 + (f2 - f1)/3;
    cutCell_.calcSubCell(celli, f3, normal);
    a3 = cutCell_.VolumeOfFluid();

    scalar f4 = f1 + 2*(f2 - f1)/3;
    cutCell_.calcSubCell(celli, f4, normal);
    scalar a4 = cutCell_.VolumeOfFluid();

    // Building and solving Vandermonde matrix equation
    scalarField a(4), f(4), C(4), fOld(4);
    {
        a[0] = a1, a[1] = a3, a[2] = a4, a[3] = a2;
        fOld[0] = f1, fOld[1] = f3, fOld[2] = f4, fOld[3] = f2;
        f[0] = 0, f[1] = (f3-f1)/(f2-f1), f[2] = (f4-f1)/(f2-f1), f[3] = 1;
        scalarSquareMatrix M(4);
        forAll(f, i)
        {
            forAll(f, j)
            {
                M[i][j] = pow(f[i], 3 - j);
            }
        }
        // C holds the 4 polynomial coefficients
        C = a;
        LUsolve(M, C);
    }

    // Finding root with Newton method

    f3 = f[1]; a3 = a[1];
    label nIter = 0;
    scalar res = mag(a3 - alpha1);
    while (res > tol && nIter < 10*maxIter)
    {
        f3 -= (C[0]*pow3(f3) + C[1]*sqr(f3) + C[2]*f3 + C[3] - alpha1)
            /(3*C[0]*sqr(f3) + 2*C[1]*f3 + C[2]);
        a3 = C[0]*pow3(f3) + C[1]*sqr(f3) + C[2]*f3 + C[3];
        res = mag(a3 - alpha1);
        nIter++;
    }
    // Scaling back to original range
    f3 = f3*(f2 - f1) + f1;

    //Check result
    label status = cutCell_.calcSubCell(celli, f3, normal);
    const scalar VOF = cutCell_.VolumeOfFluid();
    res = mag(VOF - alpha1);

    if (res > tol)
    {
    }
    else
    {
        return status;
    }

    // If tolerance not met use the secant method  with f3 as a hopefully very
    // good initial guess to crank res the last piece down below tol

    scalar x2 = f3;
    scalar g2 = VOF - alpha1;
    scalar x1 = max(1e-3*(f2 - f1), 100*SMALL);
    x1 = max(x1, f1);
    x1 = min(x1, f2);
    cutCell_.calcSubCell(celli, x1,normal);
    scalar g1 = cutCell_.VolumeOfFluid() - alpha1;

    nIter = 0;
    scalar g0(0), x0(0);
    while (res > tol && nIter < maxIter && g1 != g2)
    {
        x0 = (x2*g1 - x1*g2)/(g1 - g2);
        status = cutCell_.calcSubCell(celli, x0, normal);
        g0 = cutCell_.VolumeOfFluid() - alpha1;
        res = mag(g0);
        x2 = x1; g2 = g1;
        x1 = x0; g1 = g0;
        nIter++;
    }

    return status;
}


// ************************************************************************* //
