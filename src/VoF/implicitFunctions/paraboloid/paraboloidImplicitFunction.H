/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2019 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
                            | Copyright (C) 2019 DLR
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
    Foam::implicitFunctions::paraboloidImplicitFunction

Description
    creates a paraboloid

    Original code supplied by Henning Scheufler, DLR (2019)

SourceFiles
    paraboloidImplicitFunction.C

\*---------------------------------------------------------------------------*/

#ifndef implicitFunction_paraboloidImplicitFunction_H
#define implicitFunction_paraboloidImplicitFunction_H

#include "implicitFunction.H"
#include "mathematicalConstants.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace implicitFunctions
{

/*---------------------------------------------------------------------------*\
                 class paraboloidImplicitFunction Declaration
\*---------------------------------------------------------------------------*/

class paraboloidImplicitFunction
:
    public implicitFunction
{
private:

    // Private Member Data

        //- Coefficients of ax^2 + bx*y + cy^2
        vector coeffs_;


    // Private Member Functions

        //- No copy construct
        paraboloidImplicitFunction(const paraboloidImplicitFunction&);

        //- No bitwise assignment
        void operator=(const paraboloidImplicitFunction&);


public:

    //- Runtime type information
    TypeName("paraboloid");


    // Constructors

        //- Construct from components
        paraboloidImplicitFunction(const vector& coeffs);

        //- Construct from dictionary
        paraboloidImplicitFunction(const dictionary& dict);


    //- Destructor
    virtual ~paraboloidImplicitFunction() = default;

    virtual scalar value(const vector& p) const
    {
        return
            coeffs_.x()*sqr(p.x())
          + coeffs_.y()*p.x()*p.y()
          + coeffs_.z()*sqr(p.y())
          - p.z();
    }

    virtual vector grad(const vector& p) const
    {
        NotImplemented;
        return vector::zero;
    }

    virtual scalar distanceToSurfaces(const vector& p) const
    {
        NotImplemented;
        return 0;
    }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace implicitFunctions
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
