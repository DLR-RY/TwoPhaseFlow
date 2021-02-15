/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2019-2019 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
                            | Copyright (C) 2019-2019 DLR
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
    Foam::polyDegree1

Description
    first degree polynominal function

    c0 + c1*x + c2*y + c3*z 

    Original code supplied by Henning Scheufler, DLR (2019)

SourceFiles
    polyDegree1.C

\*---------------------------------------------------------------------------*/

#ifndef polyDegree1_H
#define polyDegree1_H


#include "FixedList.H"
#include "vector.H"
#include "vector2D.H"
#include "scalarField.H"

#include "multiDimPolyFunctions.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class polyDegree1 Declaration
\*---------------------------------------------------------------------------*/

class polyDegree1
:
    public multiDimPolyFunctions
{



private:
    // Private data


public:


    //- Runtime type information
    TypeName("polyDegree1");
    // Constructors

    //- Construct from nTerms
    polyDegree1(const Vector<label> dirs);

    //virtual ~polyDegree1();

    virtual scalar value(const vector& vec);

    virtual const scalarField& termValues(const vector& vec);

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
