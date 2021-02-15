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
    Foam::paraboloid

Description
    fitparaboloid in form f(x,y) =

    c0*x + c1*y + c2*x^2 + c3*y^2 + c4*x*y

    geomDir has to be (1 1 -1) or (1 -1 -1)
Author
    Henning Scheufler, DLR, all rights reserved.

SourceFiles
    paraboloid.C

\*---------------------------------------------------------------------------*/

#ifndef paraboloid_H
#define paraboloid_H


#include "FixedList.H"
#include "vector.H"
#include "vector2D.H"
#include "scalarField.H"

#include "multiDimPolyFunctions.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class paraboloid Declaration
\*---------------------------------------------------------------------------*/

class paraboloid
:
    public multiDimPolyFunctions
{



private:
    // Private data

    // number of polynominal term
    label nDims_;


public:

    // Constructors
    //- Runtime type information
    TypeName("paraboloid");

    //- Construct from nTerms
    paraboloid(const Vector<label> dirs);

    //virtual ~paraboloid();

    virtual scalar value(const vector& vec);

    virtual const scalarField& termValues(const vector& vec);

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
