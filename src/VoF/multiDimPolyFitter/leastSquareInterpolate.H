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
    Foam::leastSquareInterpolate

Description
    Least square interpolation scheme

Author
    Henning Scheufler, DLR, all rights reserved.

SourceFiles
    leastSquareInterpolate.C

\*---------------------------------------------------------------------------*/

#ifndef leastSquareInterpolate_H
#define leastSquareInterpolate_H

#include "fvMesh.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "OFstream.H"



#include "multiDimPolyFitter.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class leastSquareInterpolate Declaration
\*---------------------------------------------------------------------------*/
template<class T>
class leastSquareInterpolate
{



private:
    // Private data
    multiDimPolyFitter<T> polyFitter_;

    Vector<label> geomDir_;

    label nDims_;


public:

    // Constructors

        //- Construct from fvMesh and globalIndex
        leastSquareInterpolate(word functionName,const Vector<label> geomDir);


    // Member functions

        //- returns the gradient of the cell
        T interpolate
        (
            const List<vector>& positions,
            const List<T>& listValue
        );

        Map< T > interpolate
        (
            const Map <List<vector> >& positions,
            const Map <List<T> >& listValue
        );




};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
