/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2017 OpenCFD Ltd.
    Copyright (C) 2016-2017 DHI
    Copyright (C) 2018-2019 Johan Roenby
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

#include "cutCell.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

int Foam::cutCell::debug = 0;


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::cutCell::cutCell(const fvMesh& mesh)
{
    // required as otherwise setAlphaFields might not work in parallel
    mesh.C();
    mesh.V();
    mesh.Cf();
    mesh.magSf();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::cutCell::calcCellData
(
    const DynamicList<point>& cutFaceCentres,
    const DynamicList<vector>& cutFaceAreas,
    vector& subCellCentre, scalar& subCellVolume
)
{
    // Clear the fields for accumulation
    subCellCentre = Zero;
    subCellVolume = Zero;

    // first estimate the approximate cell centre as the average of
    // face centres

    vector cEst = average(cutFaceCentres);

    // Contribution to subcell centre and volume from cut faces
    forAll(cutFaceCentres, facei)
    {
        // Calculate 3*face-pyramid volume
        scalar pyr3Vol = max(
            mag(cutFaceAreas[facei] & (cutFaceCentres[facei] - cEst)), VSMALL);

        // Calculate face-pyramid centre
        vector pc = 0.75 * cutFaceCentres[facei] + 0.25 * cEst;

        // Accumulate volume-weighted face-pyramid centre
        subCellCentre += pyr3Vol * pc;

        // Accumulate face-pyramid volume
        subCellVolume += pyr3Vol;
    }

    subCellCentre /= subCellVolume;
    subCellVolume /= 3; // formula of pyramid
}


void Foam::cutCell::calcGeomDataCutFace
(
    const DynamicList<DynamicList<point>>& faceEdges,
    const vector& subCellCentre,
    vector& faceArea,
    vector& faceCentre
)
{
    // Initial guess of face centre from edge points
    point fCentre{Zero};
    label nEdgePoints{0};
    for (const DynamicList<point>& edgePoints : faceEdges)
    {
        for (const point& p : edgePoints)
        {
            fCentre += p;
            nEdgePoints++;
        }
    }
    if (nEdgePoints > 0)
    {
        fCentre /= nEdgePoints;
    }

    vector sumN{Zero};
    scalar sumA{0};
    vector sumAc{Zero};

    // calculate area and centre
    forAll(faceEdges, ei)
    {
        const DynamicList<point>& edgePoints = faceEdges[ei];
        const label nPoints = edgePoints.size();
        for (label pi = 0; pi < nPoints - 1; pi++)
        {
            const point& nextPoint = edgePoints[pi + 1];

            vector c = edgePoints[pi] + nextPoint + fCentre;
            vector n =
                (nextPoint - edgePoints[pi]) ^ (fCentre - edgePoints[pi]);
            scalar a = mag(n);

            // Edges may have different orientation
            sumN += Foam::sign(n & sumN) *  n;
            sumA += a;
            sumAc += a * c;
        }
    }

    // This is to deal with zero-area faces. Mark very small faces
    // to be detected in e.g., processorPolyPatch.
    if (sumA < ROOTVSMALL)
    {
        faceCentre = fCentre;
        faceArea = Zero;
    }
    else
    {
        faceCentre = (1.0/3.0)*sumAc/sumA;
        faceArea = 0.5*sumN;
    }

    // Check faceArea direction and change if not pointing in the subcell
    if ((faceArea & (faceCentre - subCellCentre)) >= 0)
    {
        faceArea *= (-1);
    }
}


void Foam::cutCell::calcIsoFacePointsFromEdges
(
    const vector& faceArea,
    const vector& faceCentre,
    const DynamicList<DynamicList<point>>& faceEdges,
    DynamicList<point>& facePoints
)
{
    if (mag(faceArea) == 0)
    {
        facePoints.clear();
        return;
    }
    const vector zhat = normalised(faceArea);
    vector xhat = faceEdges[0][0] - faceCentre;
    xhat = (xhat - (xhat & zhat)*zhat);
    xhat.normalise();
    if (mag(xhat) == 0)
    {
        facePoints.clear();
        return;
    }
    vector yhat = normalised(zhat ^ xhat);
    if (mag(yhat) == 0)
    {
        facePoints.clear();
        return;
    }
    yhat.normalise();

    // Calculating all intersection points
    DynamicList<point> unsortedFacePoints(3 * faceEdges.size());
    DynamicList<scalar> unsortedFacePointAngles(3 * faceEdges.size());
    for (const DynamicList<point>& edgePoints : faceEdges)
    {
        for (const point& p : edgePoints)
        {
            unsortedFacePoints.append(p);
            unsortedFacePointAngles.append
            (
                Foam::atan2
                (
                    ((p - faceCentre) & yhat),
                    ((p - faceCentre) & xhat)
                )
            );
        }
    }

    // Sorting face points by angle and inserting into facePoints
    labelList order(unsortedFacePointAngles.size());
    Foam::sortedOrder(unsortedFacePointAngles, order);
    facePoints.append(unsortedFacePoints[order[0]]);
    for (label pi = 1; pi < order.size(); ++pi)
    {
        if
        (
            mag
            (
                unsortedFacePointAngles[order[pi]]
              - unsortedFacePointAngles[order[pi - 1]]
            ) > 1e-8)
        {
            facePoints.append(unsortedFacePoints[order[pi]]);
        }
    }
}


// ************************************************************************* //
