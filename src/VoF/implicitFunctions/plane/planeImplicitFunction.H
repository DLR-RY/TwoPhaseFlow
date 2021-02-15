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
    Foam::implicitFunctions::planeImplicitFunction

Description
    creates a plane

    Original code supplied by Henning Scheufler, DLR (2019)

SourceFiles
    planeImplicitFunction.C

\*---------------------------------------------------------------------------*/

#ifndef implicitFunction_planeImplicitFunction_H
#define implicitFunction_planeImplicitFunction_H

#include "implicitFunction.H"
#include "mathematicalConstants.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace implicitFunctions
{

/*---------------------------------------------------------------------------*\
                    Class planeImplicitFunction Declaration
\*---------------------------------------------------------------------------*/

class planeImplicitFunction
:
    public implicitFunction
{
private:

    // Private Member Data

        //- ax^2 + bx*y + cy^2
        //- Reference point
        vector origin_;

        //- Normal vector
        vector normal_;


    // Private Member Functions

        //- No copy construct
        planeImplicitFunction(const planeImplicitFunction&);

        //- No assignment
        void operator=(const planeImplicitFunction&);


public:

    //- Runtime type information
    TypeName("plane");


    // Constructors

        //- Construct from components
        planeImplicitFunction(const vector& origin, const vector& normal);

        //- Construct from dictionary (used by implicitFunctions)
        planeImplicitFunction(const dictionary& dict);


    //- Destructor
    virtual ~planeImplicitFunction() = default;


    virtual scalar value(const vector& p) const
    {
        return -normal_ & (origin_ - p);
    }

    virtual vector grad(const vector& p) const
    {
        // normal_ has the length of one
        return normal_;
    }

    virtual scalar distanceToSurfaces(const vector& p) const
    {
        // normal_ has the length of one
        return mag((p - origin_) & -normal_);
    }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace implicitFunctions
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
