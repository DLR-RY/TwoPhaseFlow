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
    Foam::surfaceForces

Description
    computes the forces accting at the interace e.g. surface tension and
    gravity

SourceFiles
    surfaceForces.C

\*---------------------------------------------------------------------------*/

#ifndef surfaceForces_H
#define surfaceForces_H

#include "typeInfo.H"
#include "runTimeSelectionTables.H"
#include "volFields.H"
#include "dimensionedScalar.H"
#include "autoPtr.H"
#include "Pair.H"

#include "IOdictionary.H"
#include "volFields.H"
#include "surfaceFields.H"

#include "surfaceTensionModel.H"
#include "surfaceTensionForceModel.H"
#include "deltaFunctionModel.H"
#include "accelerationForceModel.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class surfaceForces Declaration
\*---------------------------------------------------------------------------*/

class surfaceForces
{
private:

    //-
        dictionary surfaceForcesCoeffs_;

        //- reference to volume fraction field
        const volScalarField& alpha1_;

        const fvMesh& mesh_;


        autoPtr<surfaceTensionForceModel> surfTenForceModel_;

        autoPtr<accelerationForceModel> accModel_;


    // Private Member Functions

        //- Disallow copy construct
        surfaceForces(const surfaceForces&);

        //- Disallow default bitwise assignment
        void operator=(const surfaceForces&);

public:


    // Constructors

        //- Construct from components
        surfaceForces
        (
            const volScalarField& alpha1,
            const surfaceScalarField& phi,
            const volVectorField& U,
            const dictionary& dict
        );


    //- Destructor
    virtual ~surfaceForces()
    {}


    // Member Functions
        virtual scalar capillaryDt
        (
            const volScalarField& rho1,
            const volScalarField& rho2
        );

        virtual scalar capillaryDt
        (
            const dimensionedScalar rho1,
            const dimensionedScalar rho2
        );

        dictionary& modelDict();

        const dictionary& modelDict() const;

        void correct()
        {
            surfTenForceModel_->correct();
            accModel_->correct();
        }

        const volScalarField& acc() const
        {
            return accModel_->acc();
        }

        const surfaceScalarField& accf() const
        {
            return accModel_->accf();
        }

        tmp<volScalarField> sigma() const
        {
            return surfTenForceModel_->sigma();
        }

        const surfaceScalarField& nHatf() const
        {
            return surfTenForceModel_->nHatf();
        }

        tmp<surfaceScalarField> surfaceTensionForce();

        tmp<surfaceScalarField> accelerationForce();


};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
