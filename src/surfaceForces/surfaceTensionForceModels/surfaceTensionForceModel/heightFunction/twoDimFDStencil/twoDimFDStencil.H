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
    Foam::twoDimFDStencil

Description
    adds the volume fraction field in a given direction and stores the
    height column. With these column heights teh curvature can be computed
    with finite difference operators

Author
    Henning Scheufler, DLR, all rights reserved.

SourceFiles
    twoDimFDStencil.C

\*---------------------------------------------------------------------------*/

#ifndef twoDimFDStencil_H
#define twoDimFDStencil_H

#include "fvMesh.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "OFstream.H"
#include "cutCellIso.H"
#include "scalarMatrices.H"
#include "globalIndex.H"

#include "simpleMatrix.H"
#include "HFStencil.H"

#include "contiguous.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators
class twoDimFDStencil;
class Istream;
class Ostream;

Ostream& operator<<(Ostream&, const twoDimFDStencil&);

Istream& operator>>(Istream&, twoDimFDStencil&);


/*---------------------------------------------------------------------------*\
                        Class twoDimFDStencil Declaration
\*---------------------------------------------------------------------------*/


class twoDimFDStencil
{
private:
// Private data

    //- is needed for parallel Communiction
    //  stores the current iteration the averageColValue and the current global
    //  index
    struct iterStatus
    {
        scalar avgColVal;
        label iterI;
        label gblIdx;

        iterStatus()
        :
            avgColVal(0.5),
            iterI(0),
            gblIdx(-1)
        {

        }
    };

    //- check if problem is two dimensional
    bool twoDim_;

    //- save the addedHeights;
    scalarField heights_;

    //- direction of the height x,y or z
    label dir_;

    //- global cell number of the height needed for parallel communication
    label gblcelli_;


public:

    // iterStatus one for the pos and neg direction
    Pair<iterStatus> status;

    // Constructors
        twoDimFDStencil();

        //- Construct from fvMesh and globalIndex
        twoDimFDStencil(const bool twoDim,const label dir,const label gblcelli);

        //- Construct from components
        twoDimFDStencil
        (
            const bool twoDim,
            const scalarField& heights,
            const label dir,
            const label celli
        );

        //- reference to globalNumbering
        inline label arrayToList(const label& i,const label& j);

        //- reference to globalNumbering
        inline label arrayToList3D(const label& i,const label& j,const label& k);

        //- add addColumnHeights and return next Cell
        scalar addColumnHeight(const scalarList& stencilValues);

        //- get next cells in i,j,k ordering of the local stencil
        Vector<label> nextCell(HFStencil::orientation orientation);

        //- resetColumns
        void resetColumns()
        {
            heights_ = 0;
        }

        //- get access to column heights
        scalarField& heights()
        {
            return heights_;
        }

        //- height was found
        bool foundHeight()
        {
            return mag(status[0].avgColVal - status[1].avgColVal) > 1 - 1e-6;
        }

        //- gbl number of the starting cell
        const label& gblcelli()
        {
            return gblcelli_;
        }

        //- direction of the height
        const label& direction()
        {
            return dir_;
        }

    // Member functions

        //- computes the curvature
        scalar calcCurvature(const scalar DeltaX);


        friend Ostream& operator<<(Ostream&, const twoDimFDStencil&);

        friend Istream& operator>>(Istream&, twoDimFDStencil&);

        bool operator==
        (
            const twoDimFDStencil& FDS2
        ) const;

        bool operator!=
        (
            const twoDimFDStencil& FDS2
        ) const
        {
            return !(this->operator==(FDS2));
        }

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
