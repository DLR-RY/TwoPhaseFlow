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
    Foam::implicitFunction

Description
    Base class for implicit functions

    Original code supplied by Henning Scheufler, DLR (2019)

SourceFiles
    implicitFunction.C

\*---------------------------------------------------------------------------*/

#ifndef implicitFunction_H
#define implicitFunction_H

#include "autoPtr.H"
#include "dictionary.H"
#include "vector.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class implicitFunction Declaration
\*---------------------------------------------------------------------------*/

class implicitFunction
{
    // Private Member Functions

        //- No copy construct
        implicitFunction(const implicitFunction&) = delete;

        //- No copy assignment
        void operator=(const implicitFunction&) = delete;


public:

    //- Runtime type information
    TypeName("implicitFunction");

    //- Declare run-time constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        implicitFunction,
        dict,
        (
            const dictionary& dict
        ),
        (dict)
    );

    //- Constructor
    implicitFunction();

    //- Return a reference to the selected implicitFunction
    static autoPtr<implicitFunction> New
    (
        const word& implicitFunctionType,
        const dictionary& dict
    );

    //- Destructor
    virtual ~implicitFunction() = default;


    // Member Functions

        virtual scalar value(const vector& p) const
        {
            return GREAT;
        }

        virtual vector grad(const vector& p) const
        {
            return vector::max;
        }

        virtual scalar distanceToSurfaces(const vector& p) const
        {
            return GREAT;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
