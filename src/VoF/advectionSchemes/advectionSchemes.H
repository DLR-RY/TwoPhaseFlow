/*---------------------------------------------------------------------------*\
    Modified work | Copyright (c) 2017-2019, German Aerospace Center (DLR)
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
    Foam::advectionSchemes

Description
    Base class for advection schemes

SourceFiles
    advectionSchemes.C
    newAdvectionSchemes.C

\*---------------------------------------------------------------------------*/

#ifndef advectionSchemes_H
#define advectionSchemes_H

#include "typeInfo.H"
#include "fvc.H"
#include "runTimeSelectionTables.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "dimensionedScalar.H"
#include "autoPtr.H"
#include "IOdictionary.H"

#include "reconstructionSchemes.H"



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class advectionSchemes Declaration
\*---------------------------------------------------------------------------*/

class advectionSchemes
    : public IOdictionary
{

private:

    //-
        dictionary advectionSchemesCoeffs_;

protected:

        // reference to the VoF Field
        volScalarField& alpha1_;

        // reference to the face fluxes
        const surfaceScalarField& phi_;

        // reference to the velocity field
        const volVectorField& U_;

        //- Face volumetric transport
        surfaceScalarField alphaPhi_;

        //- Pointer to reconstruction scheme
        autoPtr<reconstructionSchemes> surf_;

        // stores excution times
        vector runTime_;




    // Private Member Functions

        //- Disallow copy construct
        advectionSchemes(const advectionSchemes&);

        //- Disallow default bitwise assignment
        void operator=(const advectionSchemes&);


public:

    //-
    dictionary& modelDict();

    //-
    const dictionary& modelDict() const;


public:

    //- Runtime type information
    TypeName("advectionSchemes");


    //- Declare run-time constructor selection table
        declareRunTimeSelectionTable
        (
            autoPtr,
            advectionSchemes,
            components,
            (
                volScalarField& alpha1,
                const surfaceScalarField& phi,
                const volVectorField& U
            ),
            (alpha1, phi,U)
        );


    // Selectors

        //- Return a reference to the selected phaseChange model
        static autoPtr<advectionSchemes> New
        (
                volScalarField& alpha1,
                const surfaceScalarField& phi,
                const volVectorField& U
        );


    // Constructors

        //- Construct from components
        advectionSchemes
        (
            const word& type,
            volScalarField& alpha1,
            const surfaceScalarField& phi,
            const volVectorField& U
        );


    //- Destructor
    virtual ~advectionSchemes();


    // Member Functions

        //- advection of the interface
        virtual void advect() = 0;

        //- advection of the interface with source terms
        virtual void advect
        (
            const volScalarField::Internal& Sp,
            const volScalarField::Internal& Su
        ) = 0;

        //- reference to interface normals
        const volVectorField& normal()
        {
            return surf_->normal();
        }

        //- reference to interface centres
        const volVectorField& centre()
        {
            return surf_->centre();
        }

        reconstructionSchemes& surf()
        {
            return surf_();
        }

        //- Return mass flux
        tmp<surfaceScalarField> getRhoPhi
        (
            const dimensionedScalar rho1,
            const dimensionedScalar rho2
        ) const
        {
            return tmp<surfaceScalarField>
            (
                new surfaceScalarField
                (
                    "rhoPhi",
                    (rho1 - rho2)*alphaPhi_ + rho2*phi_
                )
            );
        }

        // Return mass flux
        tmp<surfaceScalarField> getRhoPhi
        (
            const volScalarField& rho1,
            const volScalarField& rho2
        )
        {
            return tmp<surfaceScalarField>
            (
               new surfaceScalarField
               (
                   "rhoPhi",
                   fvc::interpolate(rho1 - rho2)*alphaPhi_ + fvc::interpolate(rho2)*phi_
               )
            );
        }

        const surfaceScalarField& alphaPhi() const
        {
            return alphaPhi_;
        }

        const vector& runTime() const
        {
            return runTime_;
        }


};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
