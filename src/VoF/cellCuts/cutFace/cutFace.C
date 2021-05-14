/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
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

#include "cutFace.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

int Foam::cutFace::debug = 0;

// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * *

void Foam::cutFace::calcSubFace
(
    const label faceI,
    const scalarList& pointStatus,
    label firstFullySubmergedPoint,
    DynamicList<point>& subFacePoints,
    DynamicList<point>& surfacePoints,
    label& faceStatus,
    vector& subFaceCentre,
    vector& subFaceArea
)
{
    const pointField& points = mesh_.points();
    const face& f = mesh_.faces()[faceI];

    if (firstFullySubmergedPoint == -1) // is in gasPhase
    {
        faceStatus = 1;
        subFaceCentre = Zero;
        subFaceArea = Zero;
        return;
    }

    // loop face and append the cuts
    // loop starts at firstFullySubmergedPoint
    for
    (
        label i = firstFullySubmergedPoint;
        i < firstFullySubmergedPoint + f.size();
        ++i
    )
    {
        // max two points are appended during one cycle
        label idx = i % f.size();
        label nextIdx = (i + 1) % f.size();

        if (pointStatus[idx] > 0) // append fluid point
        {
            subFacePoints.append(points[f[idx]]);
        }
        else if (pointStatus[idx] == 0) // append cut point
        {
            subFacePoints.append(points[f[idx]]);
            surfacePoints.append(points[f[idx]]);
        }

        if
        (
            (pointStatus[idx] < 0 && pointStatus[nextIdx] > 0) ||
            (pointStatus[idx] > 0 && pointStatus[nextIdx] < 0)
        ) // cut on edge cut Value is zero
        {
            label nextP = f.nextLabel(idx);
            vector dir = points[nextP] - points[f[idx]];
            scalar weight =
                (0.0 - pointStatus[idx]) /
                (pointStatus[nextIdx] - pointStatus[idx]); // cutValue is zero

            point p = points[f[idx]] + weight * dir;

            subFacePoints.append(p);
            surfacePoints.append(p);
        }
    }

    if (subFacePoints.size() >= 3)
    {
        faceStatus = 0;
        calcSubFaceCentreAndArea(subFacePoints, subFaceCentre, subFaceArea);
    }
    else
    {
        faceStatus = -1;
    }
}


void Foam::cutFace::calcSubFace
(
    const label faceI,
    const scalarList& pointStatus,
    const scalarList& weights,
    label firstFullySubmergedPoint,
    DynamicList<point>& subFacePoints,
    DynamicList<point>& surfacePoints,
    label& faceStatus,
    vector& subFaceCentre,
    vector& subFaceArea
)
{
    const pointField& points = mesh_.points();
    const face& f = mesh_.faces()[faceI];

    if (firstFullySubmergedPoint == -1) // is in gasPhase
    {
        faceStatus = 1;
        subFaceCentre = Zero;
        subFaceArea = Zero;
        return;
    }

    // loop face and append the cuts
    // loop starts at firstFullySubmergedPoint
    for
    (
        label i = firstFullySubmergedPoint;
        i < firstFullySubmergedPoint + f.size();
        ++i
    )
    {
        // max two points are appended during one cycle
        label idx = i % f.size();
        label nextIdx = (i + 1) % f.size();

        if (pointStatus[idx] > 0) // append fluid point
        {
            subFacePoints.append(points[f[idx]]);
        }
        else if (pointStatus[idx] == 0) // append cut point
        {
            subFacePoints.append(points[f[idx]]);
            surfacePoints.append(points[f[idx]]);
        }

        if
        (
            (pointStatus[idx] < 0 && pointStatus[nextIdx] > 0) ||
            (pointStatus[idx] > 0 && pointStatus[nextIdx] < 0)
        ) // cut on edge cut Value is zero
        {
            label nextP = f.nextLabel(idx);
            vector dir = points[nextP] - points[f[idx]];

            point p = points[f[idx]] + weights[idx] * dir;

            subFacePoints.append(p);
            surfacePoints.append(p);
        }
    }

    if (subFacePoints.size() >= 3)
    {
        faceStatus = 0;
        calcSubFaceCentreAndArea(subFacePoints, subFaceCentre, subFaceArea);
    }
    else
    {
        faceStatus = -1;
    }
}


void Foam::cutFace::calcSubFace
(
    const face& f,
    const pointField& points,
    const scalarList& pointStatus,
    label firstFullySubmergedPoint,
    DynamicList<point>& subFacePoints,
    DynamicList<point>& surfacePoints,
    label& faceStatus,
    vector& subFaceCentre,
    vector& subFaceArea
)
{
    if (firstFullySubmergedPoint == -1) // in Gas
    {
        faceStatus = 1;
        subFaceCentre = Zero;
        subFaceArea = Zero;
        return;
    }

    // loop face and append the cuts
    for
    (
        label i = firstFullySubmergedPoint;
        i < firstFullySubmergedPoint + f.size();
        ++i
    )
    {
        // max two points are appended during one cycle
        label idx = i % f.size();
        label nextIdx = (i + 1) % f.size();

        if (pointStatus[idx] > 0) // append fluid point
        {
            subFacePoints.append(points[f[idx]]);
        }
        else if (pointStatus[idx] == 0) // append cut point
        {
            subFacePoints.append(points[f[idx]]);
            surfacePoints.append(points[f[idx]]);
        }

        if
        (
            (pointStatus[idx] < 0 && pointStatus[nextIdx] > 0) ||
            (pointStatus[idx] > 0 &&  pointStatus[nextIdx] < 0)
        )
        {
            label nextP = f.nextLabel(idx);
            vector dir = points[nextP] - points[f[idx]];
            scalar weight =
                (0.0 - pointStatus[idx]) /
                (pointStatus[nextIdx] - pointStatus[idx]);

            point p = points[f[idx]] + weight * dir;

            subFacePoints.append(p);
            surfacePoints.append(p);
        }
    }

    if (subFacePoints.size() >= 3)
    {
        faceStatus = 0;
        calcSubFaceCentreAndArea(subFacePoints, subFaceCentre, subFaceArea);
    }
    else
    {
        faceStatus = -1;
    }
}


void Foam::cutFace::calcSubFaceCentreAndArea
(
    DynamicList<point>& subFacePoints,
    vector& subFaceCentre,
    vector& subFaceArea
)
{
    const label nPoints = subFacePoints.size();

    // If the face is a triangle, do a direct calculation for efficiency
    // and to avoid round-off error-related problems
    if (nPoints == 3)
    {
        subFaceCentre =
                (1.0 / 3.0) * (subFacePoints[0] + subFacePoints[1] + subFacePoints[2]);

        subFaceArea = 0.5 * ((subFacePoints[1] - subFacePoints[0]) ^
                             (subFacePoints[2] - subFacePoints[0]));
    }
    else if (nPoints > 0)
    {
        vector sumN{Zero};
        scalar sumA{0};
        vector sumAc{Zero};

        point fCentre = subFacePoints[0];
        // initial guess of centre as average of subFacePoints
        for (label pi = 1; pi < nPoints; pi++)
        {
            fCentre += subFacePoints[pi];
        }

        fCentre /= nPoints;

        // loop sub triangles
        for (label pi = 0; pi < nPoints; pi++)
        {
            const point& nextPoint = subFacePoints[(pi + 1) % nPoints];

            vector c = subFacePoints[pi] + nextPoint + fCentre;
            vector n =
                (nextPoint - subFacePoints[pi]) ^ (fCentre - subFacePoints[pi]);
            scalar a = mag(n);

            sumN += n;
            sumA += a;
            sumAc += a * c;
        }

        // This is to deal with zero-area faces. Mark very small faces
        // to be detected in e.g., processorPolyPatch.
        if (sumA < ROOTVSMALL)
        {
            subFaceCentre = fCentre;
            subFaceArea = Zero;
        }
        else
        {
            subFaceCentre = (1.0 / 3.0) * sumAc / sumA;
            subFaceArea = 0.5 * sumN;
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::cutFace::cutFace(const fvMesh& mesh)
:
    mesh_(mesh)
{}


// ************************************************************************* //
