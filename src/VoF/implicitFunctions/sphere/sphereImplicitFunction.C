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

\*---------------------------------------------------------------------------*/

#include "sphereImplicitFunction.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace implicitFunctions
    {
        defineTypeNameAndDebug(sphereImplicitFunction, 0);
        addToRunTimeSelectionTable
        (
            implicitFunction,
            sphereImplicitFunction,
            dict
        );
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::implicitFunctions::sphereImplicitFunction::sphereImplicitFunction
(
    const point& origin,
    const scalar radius,
    const scalar scale
)
:
    origin_(origin),
    radius_(radius),
    scale_(scale)
{}


Foam::implicitFunctions::sphereImplicitFunction::sphereImplicitFunction
(
    const dictionary& dict
)
:
    origin_(dict.get<point>("origin")),
    radius_(dict.get<scalar>("radius")),
    scale_(dict.lookupOrDefault<scalar>("scale" ,1))
{}


// ************************************************************************* //
