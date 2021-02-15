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
    Foam::singleComponentSatProp

Description
    base class for saturation propierties

SourceFiles
    singleComponentSatProp.C
    newsingleComponentSatProp.C

\*---------------------------------------------------------------------------*/

#ifndef singleComponentSatProp_H
#define singleComponentSatProp_H

#include "typeInfo.H"
#include "runTimeSelectionTables.H"
#include "volFields.H"
#include "dimensionedScalar.H"
#include "autoPtr.H"
#include "Pair.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class singleComponentSatProp Declaration
\*---------------------------------------------------------------------------*/

class singleComponentSatProp
    : public IOobject
{

private:

    //-
        dictionary singleComponentSatPropCoeffs_;

    //-
        const fvMesh& mesh_;

    //-
        volScalarField TSat_;

    //-
        volScalarField pSat_;

    //-
        volScalarField L_;

protected:

    // Protected data

        //-
        scalar Tmin_;

        //-
        scalar Tmax_;

        //-
        dictionary& modelDict();

        //- modelDict
        const dictionary& modelDict() const;

private:

    // Private Member Functions

        //- Disallow copy construct
        singleComponentSatProp(const singleComponentSatProp&);

        //- Disallow default bitwise assignment
        void operator=(const singleComponentSatProp&);

public:

    //- Runtime type information
    TypeName("singleComponentSatProp");

    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            singleComponentSatProp,
            components,
            (
                const fvMesh& mesh,
                const dictionary& dict

            ),
            (mesh,dict)
        );

    // Selectors

        //- Return a reference to the selected phaseChange model
        static autoPtr<singleComponentSatProp> New
        (
            const fvMesh& mesh,
            const dictionary& dict
        );

    // Constructors

        //- Construct from components
        singleComponentSatProp
        (
            const fvMesh& mesh,
            const dictionary& dict
        );

    //- Destructor
    virtual ~singleComponentSatProp()
    {}

        void correct
        (
            const volScalarField& p,
            const volScalarField& T
        );
        //-
        virtual scalar satP (scalar T) const;

        //-
        virtual scalar satT (scalar p) const;

        //
        virtual scalar satL (scalar p) const;

        // //- return saturation pressure
        const volScalarField& pSat()
        {
            return pSat_;
        }

        //- return saturation Temperature
        const volScalarField& TSat()
        {
            return TSat_;
        }

        //- return latentHeat
        const volScalarField& L()
        {
            return L_;
        }

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
