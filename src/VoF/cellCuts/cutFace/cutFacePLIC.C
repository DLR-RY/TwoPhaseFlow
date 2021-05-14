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

#include "cutFacePLIC.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::cutFacePLIC::cutFacePLIC(const fvMesh& mesh)
:
    cutFace(mesh),
    mesh_(mesh),
    subFaceCentre_(Zero),
    subFaceArea_(Zero),
    subFacePoints_(10),
    surfacePoints_(4),
    pointStatus_(10),
    weight_(10),
    faceStatus_(-1)
{
    clearStorage();
}


// * * * * * * * * * * * Public Member Functions  * * * * * * * * * * * * * //

Foam::label Foam::cutFacePLIC::calcSubFace
(
    const label faceI,
    const vector& normal,
    const vector& base
 )
{
    clearStorage();

    const face& f = mesh_.faces()[faceI];
    label inLiquid = 0;
    label firstFullySubmergedPoint = -1;

    // loop face
    forAll(f, i)
    {
        // pointStatus is the distance to the plane
        scalar value = (mesh_.points()[f[i]] - base) & normal;
        if (mag(value) < SMALL)
        {
            value = 0;
        }

        pointStatus_.append(value);
        if (pointStatus_[i] > 0)
        {
            inLiquid++;
            if (firstFullySubmergedPoint == -1)
            {
                firstFullySubmergedPoint = i;
            }
        }
    }

    if (inLiquid == f.size()) // fluid face
    {
        faceStatus_ = -1;
        subFaceCentre_ = mesh_.faceCentres()[faceI];
        subFaceArea_ = mesh_.faceAreas()[faceI];
        return faceStatus_;
    }
    else if (inLiquid == 0) // gas face
    {
        faceStatus_ = 1;
        subFaceCentre_ = Zero;
        subFaceArea_ = Zero;
        return faceStatus_;
    }


    cutFace::calcSubFace
    (
        faceI,
        pointStatus_,
        firstFullySubmergedPoint,
        subFacePoints_,
        surfacePoints_,
        faceStatus_,
        subFaceCentre_,
        subFaceArea_
    );

    return faceStatus_;
}


const Foam::point& Foam::cutFacePLIC::subFaceCentre() const
{
    return subFaceCentre_;
}


const Foam::vector& Foam::cutFacePLIC::subFaceArea() const
{
    return subFaceArea_;
}


const Foam::DynamicList<Foam::point>& Foam::cutFacePLIC::subFacePoints() const
{
    return subFacePoints_;
}


const Foam::DynamicList<Foam::point>& Foam::cutFacePLIC::surfacePoints() const
{
    return surfacePoints_;
}


void Foam::cutFacePLIC::clearStorage()
{
    subFaceCentre_ = Zero;
    subFaceArea_ = Zero;
    subFacePoints_.clear();
    surfacePoints_.clear();
    pointStatus_.clear();
    weight_.clear();
    faceStatus_ = -1;
}


// ************************************************************************* //
