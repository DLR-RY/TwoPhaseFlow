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
    Foam::alphaCSF

Description

SourceFiles
    alphaCSF.C

\*---------------------------------------------------------------------------*/

#ifndef alphaCSF_H
#define alphaCSF_H

#include "typeInfo.H"
#include "volFields.H"
#include "dimensionedScalar.H"
#include "autoPtr.H"

#include "deltaFunctionModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class alphaCSF Declaration
\*---------------------------------------------------------------------------*/

class alphaCSF
    : public deltaFunctionModel
{
private:

    // Private Member Functions

        //- Disallow copy construct
        alphaCSF(const alphaCSF&);

        //- Disallow default bitwise assignment
        void operator=(const alphaCSF&);

    protected:

        //- Re-calculate the interface curvature
        virtual void calculateDeltaFunction();



public:

    //- Runtime type information
    TypeName("alphaCSF");

    // Constructors

        //- Construct from components
        alphaCSF
        (
            const dictionary& dict,
            const volScalarField& alpha1
        );

        //- Destructor
        virtual ~alphaCSF()
        {}


};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
