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
    Foam::implicitFunctions::sphereImplicitFunction

Description
    Creates a sphere

    Original code supplied by Henning Scheufler, DLR (2019)

SourceFiles
    sphereImplicitFunction.C

\*---------------------------------------------------------------------------*/

#ifndef implicitFunction_sphereImplicitFunction_H
#define implicitFunction_sphereImplicitFunction_H

#include "implicitFunction.H"
#include "point.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace implicitFunctions
{

/*---------------------------------------------------------------------------*\
                   Class sphereImplicitFunction Declaration
\*---------------------------------------------------------------------------*/

class sphereImplicitFunction
:
    public implicitFunction
{
    // Private Member Data

        //- Origin point
        const point origin_;

        //- Radius
        const scalar radius_;

        const scalar scale_;


    // Private Member Functions

        //- No copy construct
        sphereImplicitFunction(const sphereImplicitFunction&);

        //- No copy assignment
        void operator=(const sphereImplicitFunction&) = delete;


public:

    //- Runtime type information
    TypeName("sphere");


    // Constructors

        //- Construct from components
        sphereImplicitFunction
        (
            const point&,
            const scalar radius,
            const scalar scale
        );

        //- Construct from dictionary
        sphereImplicitFunction(const dictionary& dict);

    //- Destructor
    virtual ~sphereImplicitFunction() = default;


    // Member Functions

        virtual scalar value(const vector& p) const
        {
            return (-mag(p - origin_) + radius_)*scale_;
        }

        virtual vector grad(const vector& p) const
        {
            return (origin_ - p)*scale_;
        }

        virtual scalar distanceToSurfaces(const vector& p) const
        {
            return mag(mag(p-origin_) - radius_)*scale_;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace implicitFunctions
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
