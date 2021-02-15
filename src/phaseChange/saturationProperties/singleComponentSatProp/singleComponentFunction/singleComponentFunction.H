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
    Foam::singleComponentFunction

Description
    uses Function1 to compute the saturation properties

SourceFiles
    singleComponentFunction.C

\*---------------------------------------------------------------------------*/

#ifndef ClausiusClapeyron_H
#define ClausiusClapeyron_H

#include "typeInfo.H"
#include "singleComponentSatProp.H"
#include "volFields.H"
#include "dimensionedScalar.H"
#include "autoPtr.H"
#include "Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class singleComponentFunction Declaration
\*---------------------------------------------------------------------------*/

class singleComponentFunction
    : public singleComponentSatProp
{

private:

    //-
    autoPtr<Function1<scalar>>   pSatFunc_;

    //-
    autoPtr<Function1<scalar>>  TSatFunc_;

    //-
    autoPtr<Function1<scalar>>  LFunc_;



    // Private Member Functions

        //- Disallow copy construct
        singleComponentFunction(const singleComponentFunction&);

        //- Disallow default bitwise assignment
        void operator=(const singleComponentFunction&);


public:

    //- Runtime type information
    TypeName("function");


    // Constructors

        //- Construct from components
        singleComponentFunction
        (
            const fvMesh& mesh,
            const dictionary& dict
        );


    //- Destructor
    virtual ~singleComponentFunction()
    {}


    // Member Functions

        //- Read the transportProperties dictionary and update
        //-
        virtual scalar satP (scalar T) const;

        //-
        virtual scalar satT (scalar p) const;

        //
        virtual scalar satL (scalar p) const;


        // virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
