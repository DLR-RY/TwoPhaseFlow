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
    Foam::energySourceTermModel

Description
    base class of the energy source terms models

SourceFiles
    energySourceTermModel.C
    newEnergySourceTermModel.C

\*---------------------------------------------------------------------------*/

#ifndef energySourceTermModel_H
#define energySourceTermModel_H

#include "typeInfo.H"
#include "runTimeSelectionTables.H"
#include "volFields.H"
#include "dimensionedScalar.H"
#include "autoPtr.H"
#include "Pair.H"

#include "phaseModel.H"
#include "fvMatricesFwd.H"
#include "fvc.H"

#include "reconstructionSchemes.H"
#include "singleComponentSatProp.H"
#include "compressibleInterPhaseTransportModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class energySourceTermModel Declaration
\*---------------------------------------------------------------------------*/

class energySourceTermModel
    : public dictionary
{
private:
    //-
        dictionary energySourceTermModelCoeffs_;

protected:

        //- Thermo-package of phase 1
        const phaseModel& phase1_;

        //- Thermo-package of phase 2
        const phaseModel& phase2_;

        //- reference to turbulence Model
        const compressibleInterPhaseTransportModel& turbModel_;

        //- reference to pressure
        const volScalarField& p_;

        //- reference to saturation model
        singleComponentSatProp& satModel_;

        //- reference to reconstruction scheme
        reconstructionSchemes& surf_;

    // Private Member Functions

        //- Disallow copy construct
        energySourceTermModel(const energySourceTermModel&);

        //- Disallow default bitwise assignment
        void operator=(const energySourceTermModel&);

public:

    //-
    dictionary& modelDict();

    //-
    const dictionary& modelDict() const;

public:

    //- Runtime type information
    TypeName("energySourceTermModel");

    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            energySourceTermModel,
            components,
            (
                const phaseModel& phase1,
                const phaseModel& phase2,
                const compressibleInterPhaseTransportModel& turbModel,
                const volScalarField& p,
                singleComponentSatProp& satModel,
                reconstructionSchemes& surf,
                const dictionary& dict
            ),
            (phase1,phase2,turbModel,p,satModel,surf,dict)
        );

    // Selectors

        //- Return a reference to the selected phaseChange model
        static autoPtr<energySourceTermModel> New
        (
            const phaseModel& phase1,
            const phaseModel& phase2,
            const compressibleInterPhaseTransportModel& turbModel,
            const volScalarField& p,
            singleComponentSatProp& satModel,
            reconstructionSchemes& surf,
            const dictionary& dict
        );

    // Constructors

        //- Construct from components
        energySourceTermModel
        (
            const word& type,
            const phaseModel& phase1,
            const phaseModel& phase2,
            const compressibleInterPhaseTransportModel& turbModel,
            const volScalarField& p,
            singleComponentSatProp& satModel,
            reconstructionSchemes& surf,
            const dictionary& dict
        );

    //- Destructor
    virtual ~energySourceTermModel()
    {}

    // Member Functions

        virtual tmp<fvScalarMatrix> TSource1() = 0;

        virtual tmp<fvScalarMatrix> TSource2() = 0;

        virtual tmp<volScalarField> energySource() = 0;

        virtual tmp<volScalarField> energyFlux1() = 0;

        virtual tmp<volScalarField> energyFlux2() = 0;

        virtual tmp<surfaceScalarField> getalphaf()
        {
            return fvc::interpolate(phase1_);
        }

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
