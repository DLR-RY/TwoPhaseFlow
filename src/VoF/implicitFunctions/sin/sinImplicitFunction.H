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
    Foam::implicitFunctions::sinImplicitFunction

Description

    Original code supplied by Henning Scheufler, DLR (2019)

SourceFiles
    sinImplicitFunction.C

\*---------------------------------------------------------------------------*/

#ifndef implicitFunction_sinImplicitFunction_H
#define implicitFunction_sinImplicitFunction_H

#include "implicitFunction.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace implicitFunctions
{

/*---------------------------------------------------------------------------*\
                     Class sinImplicitFunction Declaration
\*---------------------------------------------------------------------------*/

class sinImplicitFunction
:
    public implicitFunction
{
private:

    // Private Member Data

        //- origin point
        const scalar period_;

        //- Radius
        const scalar phase_;

        const scalar amplitude_;

        vector up_;

        vector direction_;

        vector origin_;


    // Private Member Functions

        //- No copy construct
        sinImplicitFunction(const sinImplicitFunction&) = delete;

        //- No bitwise assignment
        void operator=(const sinImplicitFunction&) = delete;


public:

    //- Runtime type information
    TypeName("sin");


    // Constructors

        //- Construct from components
        sinImplicitFunction
        (
            const scalar period,
            const scalar phase,
            const scalar amplitude,
            const vector& direction,
            const vector& up,
            const vector& origin
        );

        //- Construct from dictionary (used by implicitFunctions)
        sinImplicitFunction(const dictionary& dict);


    //- Destructor
    virtual ~sinImplicitFunction() = default;

    virtual scalar value(const vector& p) const
    {
        const scalar x = (p - origin_) & direction_;
        const scalar z = (p - origin_) & -up_;

        return
            z
          - amplitude_
           *Foam::sin(2*constant::mathematical::pi*x/period_ + phase_);
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
