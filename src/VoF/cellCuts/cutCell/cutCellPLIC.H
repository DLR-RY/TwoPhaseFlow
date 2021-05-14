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

Class
    Foam::cutCellPLIC

Description
    Class for cutting a cell, cellI, of an fvMesh, mesh_, at its intersection
    with an surface defined by a normal and cutValue_ (defined as distance to
    the cell centre).

    Reference:
    \verbatim
        Henning Scheufler, Johan Roenby,
        Accurate and efficient surface reconstruction from volume
        fraction data on general meshes,
        Journal of Computational Physics, 2019,
        doi 10.1016/j.jcp.2019.01.009

    \endverbatim


    Original code supplied by Henning Scheufler, DLR (2019)

SourceFiles
    cutCellPLIC.C

\*---------------------------------------------------------------------------*/

#ifndef cutCellPLIC_H
#define cutCellPLIC_H

#include "cutCell.H"
#include "cutFacePLIC.H"
#include "fvMesh.H"
#include "surfaceFields.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class cutCellPLIC Declaration
\*---------------------------------------------------------------------------*/

class cutCellPLIC
:
    public cutCell
{
    // Private Data

        //- Mesh whose cells and faces to cut at their intersection with an
        //  isosurface.
        const fvMesh& mesh_;

        //- Cell to cut
        label cellI_;

        //- Normal of the cutting plane
        vector normal_;

        //- Cutvalue used to cut cell
        scalar cutValue_;

        //- An cutFacePLIC object to get access to its face cutting functionality
        cutFacePLIC cutFace_;

        //- List of face centres for CutFaces
        DynamicList<point> cutFaceCentres_;

        //- List of face area vectors for PLICCutFaces
        DynamicList<vector> cutFaceAreas_;

        //- Storage for subFace edges belonging to PLICFace
        DynamicList<DynamicList<point>> plicFaceEdges_;

        //- Points constituting the cell-PLICsurface intersection (isoface)
        DynamicList<point> facePoints_;

        //- Face centre of the cutFace
        point faceCentre_;

        //- Face normal of the PLICface by convention pointing from high to low
        //  values (i.e. opposite of the gradient vector).
        vector faceArea_;

        //- Cell centre of the subcell of celli which is "fully submerged", i.e.
        //  where the function value is higher than the isoValue_
        point subCellCentre_;

        //- Volume of fully submerged subcell
        scalar subCellVolume_;

        //- Volume of Fluid for cellI (subCellVolume_/mesh_.V()[cellI])
        scalar VOF_;

        //- A cell status label taking one of the values:
        //
        //   -1: cell is fully below the PLICsurface
        //    0: cell is cut
        //   +1: cell is fully above the PLICsurface
        label cellStatus_;


  public:

    // Constructors

        //- Construct from fvMesh
        explicit cutCellPLIC(const fvMesh& mesh);


    // Member Functions

        //- Sets internal values and returns face status
        label calcSubCell
        (
            const label celli,
            const scalar cutValue,
            const vector& normal
        );

        //- Returns subCellCentre
        const point& subCellCentre() const;

        //- Returns subCellVolume
        scalar subCellVolume() const;

        //- Returns the points of the cutting PLICface
        const DynamicList<point>& facePoints();

        //- Returns the centre of the cutting PLICface
        const point& faceCentre() const;

        //- Returns the area normal vector of the cutting PLICface
        const vector& faceArea() const;

        //- Returns cellStatus
        label cellStatus() const;

        //- Returns volume of fluid value
        scalar VolumeOfFluid() const;

        //- Returns cutValue
        scalar cutValue() const;

        //- Resets internal values
        void clearStorage();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
