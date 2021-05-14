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

#include "cutCellPLIC.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::cutCellPLIC::cutCellPLIC(const fvMesh& mesh)
:
    cutCell(mesh),
    mesh_(mesh),
    cellI_(-1),
    normal_(Zero),
    cutValue_(0),
    cutFace_(mesh_),
    cutFaceCentres_(10),
    cutFaceAreas_(10),
    plicFaceEdges_(10),
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

Foam::label Foam::cutCellPLIC::calcSubCell
(
    const label celli,
    const scalar cutValue,
    const vector& normal
)
{
    clearStorage();
    normal_ = normal;
    cellI_ = celli;
    cutValue_ = cutValue;
    const cell& c = mesh_.cells()[celli];

    vector base = mesh_.C()[cellI_] + normal_ * cutValue_;
    bool fullyBelow = true;
    bool fullyAbove = true;

    label nFaceBelowInterface = 0;

    // loop over all cell faces
    for (const label faceI : c)
    {
        const label faceStatus = cutFace_.calcSubFace(faceI, normal_, base);

        if (faceStatus == 0) // face is cut
        {
            cutFaceCentres_.append(cutFace_.subFaceCentre());
            cutFaceAreas_.append(cutFace_.subFaceArea());
            plicFaceEdges_.append(cutFace_.surfacePoints());
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
            plicFaceEdges_,
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


const Foam::point& Foam::cutCellPLIC::subCellCentre() const
{
    return subCellCentre_;
}


Foam::scalar Foam::cutCellPLIC::subCellVolume() const
{
    return subCellVolume_;
}


const Foam::DynamicList<Foam::point>& Foam::cutCellPLIC::facePoints()
{
    if (facePoints_.size() == 0)
    {
        // get face points in sorted order
        calcIsoFacePointsFromEdges
        (
            faceArea_,
            faceCentre_,
            plicFaceEdges_,
            facePoints_
        );
    }

    return facePoints_;
}


const Foam::point& Foam::cutCellPLIC::faceCentre() const
{
    return faceCentre_;
}


const Foam::vector& Foam::cutCellPLIC::faceArea() const
{
    return faceArea_;
}


Foam::scalar Foam::cutCellPLIC::VolumeOfFluid() const
{
    return VOF_;
}


Foam::label Foam::cutCellPLIC::cellStatus() const
{
    return cellStatus_;
}


Foam::scalar Foam::cutCellPLIC::cutValue() const
{
    return cutValue_;
}


void Foam::cutCellPLIC::clearStorage()
{
    cellI_ = -1;
    cutValue_ = 0;
    cutFaceCentres_.clear();
    cutFaceAreas_.clear();
    plicFaceEdges_.clear();
    facePoints_.clear();
    faceCentre_ = Zero;
    faceArea_ = Zero;
    subCellCentre_ = Zero;
    subCellVolume_ = -10;
    VOF_ = -10;
    cellStatus_ = -1;
}


// ************************************************************************* //
