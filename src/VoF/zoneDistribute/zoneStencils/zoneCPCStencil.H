/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
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
    Foam::zoneCPCStencil

Description
    computes a cell point cell stencil in a narrow band. resizes in case
    of topological change

    Original code supplied by Henning Scheufler, DLR (2019)

SourceFiles
    zoneCPCStencil.C

\*---------------------------------------------------------------------------*/

#ifndef zoneCPCStencil_H
#define zoneCPCStencil_H

#include "zoneCellStencils.H"
#include "boolList.H"
#include "HashSet.H"
#include "Map.H"
#include "indirectPrimitivePatch.H"
#include "fvMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class zoneCPCStencil Declaration
\*---------------------------------------------------------------------------*/

class zoneCPCStencil
:
    public MeshObject
    <
        fvMesh,
        TopologicalMeshObject,
        zoneCPCStencil
    >,
    public zoneCellStencils
{
    // Private Data

        labelList nonEmptyBoundaryPoints_;

        //- is valid boundaryface ?
        //- non valid face are on empty, wedge and processor patches
        boolList isValidBFace_;

        //- cell is already uptodate
        boolList uptodate_;


    // Private Member Functions

        //- Sync boundary points
        Map<bool> syncCoupledBoundaryPoints
        (
            const boolList& zone,
            const labelList& boundaryPoints
        ) const;

        //- Calculates per point the neighbour data (= pointCells)
        void calcPointBoundaryData
        (
            const boolList& zone,
            const boolList& isValidBFace,
            const labelList& boundaryPoints,
            Map<labelList>& neiGlobal
        ) const;

        //- Calculate stencil
        void calculateStencil
        (
            const boolList& zone,
            labelListList& globalCellCells
        );

        //- No copy construct
        zoneCPCStencil(const zoneCPCStencil&) = delete;

        //- No copy assignment
        void operator=(const zoneCPCStencil&) = delete;

public:

    //- Runtime information
    TypeName("zoneCPCStencil");


    // Constructors

        //- Construct from all cells and boundary faces
        explicit zoneCPCStencil(const fvMesh&);

    // Selectors
        static zoneCPCStencil& New(const fvMesh&);

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
