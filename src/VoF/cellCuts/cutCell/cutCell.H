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

Class
    Foam::cutCell

Description
    Service routines for cutting a cell, celli, of an fvMesh, mesh_,
    at its intersection with a surface

    Original code supplied by
        Johan Roenby, DHI (2016)
        Henning Scheufler, DLR (2019)

SourceFiles
    cutCell.C

\*---------------------------------------------------------------------------*/

#ifndef cutCell_H
#define cutCell_H

#include "fvMesh.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class cutCell Declaration
\*---------------------------------------------------------------------------*/

class cutCell
{
protected:

    // Protected Member Functions

        //- Calculates volume and centre of the cutted cell
        static void calcCellData
        (
            const DynamicList<point>& cutFaceCentres,
            const DynamicList<vector>& cutFaceAreas,
            vector& subCellCentre,
            scalar& subCellVolume
        );

        //- Calculates area and centre of the cutting face
        static void calcGeomDataCutFace
        (
            const DynamicList<DynamicList<point>>& faceEdges,
            const vector& subCellCentre,
            vector& faceArea,
            vector& faceCentre
        );

        //- Calculates the point of the cutting face
        static void calcIsoFacePointsFromEdges
        (
            const vector& faceArea,
            const vector& faceCentre,
            const DynamicList<DynamicList<point>> & faceEdges,
            DynamicList<point>& facePoints
        );

public:

    // Static Data
    static int debug;


    // Constructors

        //- Construct from fvMesh
        explicit cutCell(const fvMesh& mesh);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
