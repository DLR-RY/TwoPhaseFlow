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
    Foam::zoneCellStencils

Description
    base class for cell stencil in a narrow band

    Original code supplied by Henning Scheufler, DLR (2019)

SourceFiles
    zoneCellStencils.C

\*---------------------------------------------------------------------------*/

#ifndef zoneCellStencils_H
#define zoneCellStencils_H

#include "boolList.H"
#include "HashSet.H"
#include "Map.H"
#include "indirectPrimitivePatch.H"
#include "MeshObject.H"
#include "fvMesh.H"
#include "globalIndex.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class zoneCellStencils Declaration
\*---------------------------------------------------------------------------*/

class zoneCellStencils
:
    public labelListList
{
protected:

    // Protected Members

        const fvMesh& meshRef_;

        labelHashSet needComm_;

        globalIndex globalNumbering_;


    // Protected Member Functions

        //- Return patch of all coupled faces.
        autoPtr<indirectPrimitivePatch> nonEmptyFacesPatch() const;

        //- Return patch of all coupled faces.
        autoPtr<indirectPrimitivePatch> allCoupledFacesPatch() const;

        //- Valid boundary faces (not empty and not coupled)
        void validBoundaryFaces(boolList& isValidBFace) const;

        //- Merge two lists and guarantee globalI is first.
        static void merge
        (
            const label globalI,
            const labelList& pGlobals,
            labelList& cCells
        );

        //- Combine operator for labelLists
        class unionEqOp
        {
            public:
            void operator()(labelList& x, const labelList& y) const;
        };

        //- Collect cell neighbours of faces in global numbering
        void insertFaceCells
        (
            const label exclude0,
            const label exclude1,
            const boolList& nonEmptyFace,
            const labelList& faceLabels,
            labelHashSet& globals
        ) const;

        //- Collect cell neighbours of faces in global numbering
        labelList calcFaceCells
        (
            const boolList& nonEmptyFace,
            const labelList& faceLabels,
            labelHashSet& globals
        ) const;

        virtual void calculateStencil
        (
            const boolList& zone,
            labelListList& globalCellCells
        ) = 0;


        //- No copy construct
        zoneCellStencils(const zoneCellStencils&) = delete;

        //- No copy assignment
        void operator=(const zoneCellStencils&) = delete;

public:

        // Declare name of the class and its debug switch
        TypeName("zoneCellStencils");


    // Constructors

        //- Construct from all cells and boundary faces
        explicit zoneCellStencils(const fvMesh&);

        //- Calculates per cell the neighbour data
        //  (= cell or boundary in global numbering).
        //  First element is always cell itself!
        //- Destructor
        virtual ~zoneCellStencils() = default;


    // Member Functions

        void updateStencil
        (
            const boolList& zone
        )
        {
            calculateStencil(zone,*this);
        }

        const labelHashSet& needsComm()
        {
            return needComm_;
        }

        const polyMesh& mesh() const
        {
            return meshRef_;
        }

         //- Global numbering for cells and boundary faces
        const globalIndex& globalNumbering() const
        {
             return globalNumbering_;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
