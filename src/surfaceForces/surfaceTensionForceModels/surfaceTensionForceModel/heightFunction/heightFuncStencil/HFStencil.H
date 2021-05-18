/*---------------------------------------------------------------------------*\
            Copyright (c) 2017-2019, German Aerospace Center (DLR)
-------------------------------------------------------------------------------
License
    This file is part of the VoFLibrary source code library, which is an
	unofficial extension to OpenFOAM.
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
    Foam::HFStencil

Description
    Class that provides cell-point-cell stencil. A labelList<List <label>>
    is returned. The elements are in globalCellNumbering (cell and boundaryFaces)
    (see. class globalIndex). The boundary values and processor information are
    provided by a Map. Assumes a hex grid resulting 27 cells in 3D and 9 Cells
    in 2D. The cells are ordered (i,j,k).

    Assumption i j k have to be alligned with x y and z:
        -> mesh rotation would break the method

Author
    Henning Scheufler, DLR, all rights reserved.

SourceFiles
    HFStencil.C

\*---------------------------------------------------------------------------*/

#ifndef HFStencil_H
#define HFStencil_H

#include "OFstream.H"
#include "fvMesh.H"
#include "globalIndex.H"
#include "token.H"
#include "volFields.H"
#include "zoneDistribute.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class HFStencil Declaration
\*---------------------------------------------------------------------------*/

class HFStencil
{
        //- Reference to mesh
        const fvMesh& mesh_;

        //- cell-point-cell stencil elements are in global addressing
        //  2D Case stencil has a size of 9
        //  3D Case stencil has a size of 27
        List<labelList > stencilHF_;

        //- check wether the cell is a cuboid or not
        boolList isCuboid_;

        //- is twodimensional
        bool twoDim_;

        //- angle tolerance is zero on a perfect grid
        scalar toleranceAngle_;

        //- sets the cuboids
        void findCuboids();

        void sortStencilToIJKFormat
        (
            const label celli,
            const labelList& CPCstencil,
            const vectorField& cellCentres,
            const vector cc
        );

        label calcPosInStencil(const vector& cc, const vector& neic);

        inline label categorizePoint(const scalar pos);

        //- reference to globalNumbering
        inline label arrayToList3D(const label i,const label j,const label k);


  public:
    // Constructors
        //- Component labeling enumeration
        enum direction { I, J, K };

        enum orientation {pos,neg};

        //- Construct from fvMesh
        HFStencil(const fvMesh&,const scalar angleTolerance);

        // Member functions

        //- updates the stencil with boolList the size has to match
        //  mesh.nCells()
        void updateStencil(const boolList& nextToInterface);

        //- return in the index of the selecet position in the stencil
        //  the values of posInStencil between 0 and 2 and for twoDim
        //  cases k or z have to be 0
        label getCellLabel
        (
            const label celli,
            const Vector<label>& posInStencil
        );

        //- return stencil reference
        const List<labelList>& getStencil()
        {
            return stencilHF_;
        }

        //- is the cell cuboid?
        const boolList& isCuboid()
        {
            return isCuboid_;
        }

        //- return addressing reference
        // const globalIndex& globalNumbering() const
        // {
        //     return exchangeFields_.globalNumbering();
        // }

        //- Returns stenicl and provides a Map with globalNumbering
        Map< scalar > getDatafromOtherProc
        (
            const boolList& nextToInterface,
            const volScalarField& alpha
        );

        // //- get non-const reference to CPC Stencil
        // zoneDistribute& getCPCStencil()
        // {
        //     return exchangeFields_;
        // }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


#endif

// ************************************************************************* //
