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
    Foam::surfaceTensionForceModel

Description

SourceFiles
    surfaceTensionForceModel.C
    newsurfaceTensionForceModel.C

\*---------------------------------------------------------------------------*/

#ifndef surfaceTensionForceModel_H
#define surfaceTensionForceModel_H

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
#include "deltaFunctionModel.H"
#include "fvc.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class surfaceTensionForceModel Declaration
\*---------------------------------------------------------------------------*/

class surfaceTensionForceModel
{

protected:

        //- reference to volume fraction field
        const volScalarField& alpha1_;

        //- reference to face velocity
        const surfaceScalarField& phi_;

        //- reference to velcoity
        const volVectorField& U_;

        autoPtr<surfaceTensionModel> surfTenModel_;

        autoPtr<deltaFunctionModel> deltaFunctionModel_;

        //- surface face normal on the face
        surfaceScalarField nHatf_;

        //- curvature field
        volScalarField K_;

        //- curvature field on faces
        surfaceScalarField Kf_;

private:

    // Private Member Functions

        //- Disallow copy construct
        surfaceTensionForceModel(const surfaceTensionForceModel&);

        //- Disallow default bitwise assignment
        void operator=(const surfaceTensionForceModel&);

public:

    //- Runtime type information
    TypeName("surfaceTensionForceModel");

    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            surfaceTensionForceModel,
            components,
            (
                const dictionary& dict,
                const volScalarField& alpha1,
                const surfaceScalarField& phi,
                const volVectorField& U
            ),
            (dict,alpha1,phi,U)
        );

    // Selectors

        //- Return a reference to the selected phaseChange model
        static autoPtr<surfaceTensionForceModel> New
        (
            const dictionary& dict,
            const volScalarField& alpha1,
            const surfaceScalarField& phi,
            const volVectorField& U
        );

    // Constructors

        //- Construct from components
        surfaceTensionForceModel
        (
            const word& type,
            const dictionary& dict,
            const volScalarField& alpha1,
            const surfaceScalarField& phi,
            const volVectorField& U
        );

    //- Destructor
    virtual ~surfaceTensionForceModel()
    {}

    // Member Functions

        const volScalarField& K() const
        {
            return K_;
        }

        const surfaceScalarField& Kf() const
        {
            return Kf_;
        }

        const surfaceScalarField& nHatf() const
        {
            return nHatf_;
        }

        tmp<volScalarField> sigma() const
        {
            return surfTenModel_->sigma();
        }

        virtual void correct() = 0;

        virtual tmp<surfaceScalarField> surfaceTensionForce()
        {
            return fvc::interpolate(sigma())*Kf()
                   * deltaFunctionModel_->deltaFunction();
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
