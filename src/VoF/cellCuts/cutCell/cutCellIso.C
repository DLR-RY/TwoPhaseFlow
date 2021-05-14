/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2017 DHI
    Copyright (C) 2018-2019 Johan Roenby
    Copyright (C) 2020 DLR
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

#include "cutCellIso.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::cutCellIso::cutCellIso(const fvMesh& mesh, scalarField& f)
:
    cutCell(mesh),
    mesh_(mesh),
    cellI_(-1),
    f_(f),
    cutValue_(0),
    cutFace_(cutFaceIso(mesh_, f_)),
    cutFaceCentres_(10),
    cutFaceAreas_(10),
    isoFaceEdges_(10),
    facePoints_(10),
    faceCentre_(Zero),
    faceArea_(Zero),
    subCellCentre_(Zero),
    subCellVolume_(-10),
    VOF_(-10),
    cellStatus_(-1)
{
    clearStorage();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::label Foam::cutCellIso::calcSubCell
(
    const label cellI,
    const scalar cutValue
)
{
    // resets data members
    clearStorage();
    cellI_ = cellI;
    cutValue_ = cutValue;
    const cell& c = mesh_.cells()[cellI_];

    bool fullyBelow = true;
    bool fullyAbove = true;

    label nFaceBelowInterface = 0;

    // loop over cell faces
    for (const label facei : c)
    {
        const label faceStatus = cutFace_.calcSubFace(facei, cutValue_);

        if (faceStatus == 0) // face is cut
        {
            cutFaceCentres_.append(cutFace_.subFaceCentre());
            cutFaceAreas_.append(cutFace_.subFaceArea());
            isoFaceEdges_.append(cutFace_.surfacePoints());
            fullyBelow = false;
            fullyAbove = false;
        }
        else if (faceStatus == -1) // face fully below
        {
            cutFaceCentres_.append(cutFace_.subFaceCentre());
            cutFaceAreas_.append(cutFace_.subFaceArea());
            fullyAbove = false;
            nFaceBelowInterface++;
        }
        else
        {
            fullyBelow = false;
        }
    }

    if (!fullyBelow && !fullyAbove) // cell cut at least at one face
    {
        cellStatus_ = 0;

        // calc faceArea and faceCentre
        calcGeomDataCutFace
        (
            isoFaceEdges_,
            average(cutFaceCentres_),
            faceArea_,
            faceCentre_
        );

        // In the rare but occuring cases where a cell is only touched at a
        // point or a line the isoFaceArea_ will have zero length and here the
        // cell should be treated as either completely empty or full.
        if (mag(faceArea_) < 10*SMALL)
        {
            if (nFaceBelowInterface == 0)
            {
                // Cell fully above isosurface
                cellStatus_ = 1;
                subCellCentre_ = Zero;
                subCellVolume_ = 0;
                VOF_ = 0;
                return cellStatus_;
            }
            else
            {
                // Cell fully below isosurface
                cellStatus_ = -1;
                subCellCentre_ = mesh_.C()[cellI_];
                subCellVolume_ = mesh_.V()[cellI_];
                VOF_ = 1;
                return cellStatus_;
            }
        }

        cutFaceCentres_.append(faceCentre_);
        cutFaceAreas_.append(faceArea_);

        // calc volume and sub cell centre
        calcCellData
        (
            cutFaceCentres_,
            cutFaceAreas_,
            subCellCentre_,
            subCellVolume_
        );

        VOF_ = subCellVolume_ / mesh_.V()[cellI_];
    }
    else if (fullyAbove) // cell fully above isosurface
    {
        cellStatus_ = 1;
        subCellCentre_ = Zero;
        subCellVolume_ = 0;
        VOF_ = 0;
    }
    else if (fullyBelow) // cell fully below isosurface
    {
        cellStatus_ = -1;
        subCellCentre_ = mesh_.C()[cellI_];
        subCellVolume_ = mesh_.V()[cellI_];
        VOF_ = 1;
    }

    return cellStatus_;
}


const Foam::point& Foam::cutCellIso::subCellCentre() const
{
    return subCellCentre_;
}


Foam::scalar Foam::cutCellIso::subCellVolume() const
{
    return subCellVolume_;
}


const Foam::DynamicList<Foam::point>& Foam::cutCellIso::facePoints()
{
    if (facePoints_.size() == 0)
    {
        // get face points in sorted order
        calcIsoFacePointsFromEdges
        (
            faceArea_,
            faceCentre_,
            isoFaceEdges_,
            facePoints_
        );
    }

    return facePoints_;
}


const Foam::point& Foam::cutCellIso::faceCentre() const
{
    return faceCentre_;
}


const Foam::vector& Foam::cutCellIso::faceArea() const
{
    return faceArea_;
}


Foam::label Foam::cutCellIso::cellStatus() const
{
    return cellStatus_;
}


Foam::scalar Foam::cutCellIso::VolumeOfFluid() const
{
    return VOF_;
}


Foam::scalar Foam::cutCellIso::cutValue() const
{
    return cutValue_;
}


void Foam::cutCellIso::clearStorage()
{
    cellI_ = -1;
    cutValue_ = 0;
    cutFaceCentres_.clear();
    cutFaceAreas_.clear();
    isoFaceEdges_.clear();
    facePoints_.clear();
    faceCentre_ = Zero;
    faceArea_ = Zero;
    subCellCentre_ = Zero;
    subCellVolume_ = -10;
    VOF_ = -10;
    cellStatus_ = -1;
}


// ************************************************************************* //
