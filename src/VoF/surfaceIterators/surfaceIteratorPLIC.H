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

Class
    Foam::surfaceIteratorPLIC

Description
    Finds the cutValue that matches the volume fraction

    Reference:
    \verbatim
        Henning Scheufler, Johan Roenby,
        Accurate and efficient surface reconstruction from volume
        fraction data on general meshes,
        Journal of Computational Physics, 2019,
        doi 10.1016/j.jcp.2019.01.009
    \endverbatim

Author
    Johan Roenby, DHI, all rights reserved.
    Modified Henning Scheufler, DLR

SourceFiles
    surfaceIteratorPLIC.C

\*---------------------------------------------------------------------------*/

#ifndef surfaceIteratorPLIC_H
#define surfaceIteratorPLIC_H

#include "fvMesh.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "cutCellPLIC.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class surfaceIteratorPLIC Declaration
\*---------------------------------------------------------------------------*/

class surfaceIteratorPLIC
{
    // Private Data

        //- Mesh whose cells and faces to cut at their intersection with an
        //- isosurface.
        const fvMesh& mesh_;

        //- Cuts the cell and returns volume centre and normal
        cutCellPLIC cutCell_;

        //- Tolerance for marking of surface cells:
        //  Those with surfCellTol_ < alpha1 < 1 - surfCellTol_
        scalar surfCellTol_;


public:

    // Constructors

        //- Construct from fvMesh and a scalarField
        //  Length of scalarField should equal number of mesh points
        surfaceIteratorPLIC(const fvMesh& mesh, const scalar tol);


    // Member Functions

        //- Determine if a cell is a surface cell
        bool isASurfaceCell(const scalar alpha1) const
        {
            return
            (
                surfCellTol_ < alpha1
             && alpha1 < 1 - surfCellTol_
            );
        }

        //- Finds matching cutValue for the given value fraction
        //  \return the cellStatus
        label vofCutCell
        (
            const label celli,
            const scalar alpha1,
            const scalar tol,
            const label maxIter,
            vector normal
        );

        //- The centre point of cutted volume
        const point& subCellCentre() const
        {
            return cutCell_.subCellCentre();
        }

        //- The volume of cutted volume
        scalar subCellVolume() const
        {
            return cutCell_.subCellVolume();
        }

        //- The centre of cutting isosurface
        const point& surfaceCentre() const
        {
            return cutCell_.faceCentre();
        }

        //- The area vector of cutting isosurface
        const vector& surfaceArea() const
        {
            return cutCell_.faceArea();
        }

        //- Volume of Fluid for cellI (subCellVolume_/mesh_.V()[cellI])
        scalar VolumeOfFluid() const
        {
            return cutCell_.VolumeOfFluid();
        }

        //- The cutValue
        scalar cutValue() const
        {
            return cutCell_.cutValue();
        }

        //- The cellStatus
        label cellStatus()
        {
            return cutCell_.cellStatus();
        }

        //- The points of the cutting isosurface in sorted order
        const DynamicList<point>& facePoints()
        {
            return cutCell_.facePoints();
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
