/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is not part of OpenFOAM.

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

Application
    generateU
    
Description
    Generates velocity field for the classical test case with a sphere 
    deformed into a shape with long tongues and back again.

Author
    Johan Roenby, DHI, all rights reserved.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    Info<< "Reading field U\n" << endl;
        
    volVectorField U
    (
        IOobject
        (
            "U",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh
    );

    scalar PI = Foam::constant::mathematical::pi;
    {
        scalarField X(mesh.C().component(0));
        scalarField Z(mesh.C().component(2));
        scalarField u(-sin(2.0*PI*Z)*pow(sin(PI*X),2.0));
        scalarField w(sin(2.0*PI*X)*pow(sin(PI*Z),2.0));
        forAll(U,ci)
        {
            U[ci] = u[ci]*vector(1.0,0.0,0.0) + 0.0*vector(0.0,1.0,0.0) + w[ci]*vector(0.0,0.0,1.0);
        }
    }
    
    Info<< "Reading/calculating face flux field phi\n" << endl;

    surfaceScalarField phi
    (
        IOobject
        (
            "phi",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        linearInterpolate(U) & mesh.Sf()
    );

    {
        const scalarField x(mesh.points().component(0));
        const scalarField y(mesh.points().component(1));
        const scalarField z(mesh.points().component(2));
        scalarField psi((1.0/PI)*pow(sin(PI*x),2.0)*pow(sin(PI*z),2.0));
        forAll(phi,fi)
        {
            phi[fi] = 0.0;
            labelList faceLabels = mesh.faces()[fi];
            label nPoints = faceLabels.size();
            forAll(faceLabels,pi)
            {
                label pl1 = faceLabels[pi];
                label pl2 = faceLabels[(pi+1) % nPoints];
                phi[fi] += 0.5*(psi[pl1]+psi[pl2])*(y[pl2]-y[pl1]);
            }
        }
    }
    scalarField sumPhi = fvc::surfaceIntegrate(phi);    
    scalar maxMagSumPhi(0.0);
    label maxLabel(0);
    forAll(sumPhi,ci)
    {
        scalar msp = mag(sumPhi[ci]);
        if (msp > maxMagSumPhi )
        {
            maxMagSumPhi = msp;
            maxLabel = ci;
        }
    }
    Info << "maxMagSumPhi/cellVol = " << maxMagSumPhi/mesh.V()[maxLabel] << endl;
    
    ISstream::defaultPrecision(18);

    U.write();
    phi.write();
    
    Info<< "End\n" << endl;

    return 0;
}

// ************************************************************************* //
