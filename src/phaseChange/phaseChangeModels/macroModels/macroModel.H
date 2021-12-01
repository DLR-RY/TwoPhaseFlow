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
    Foam::macroModel

Description
    base class of the energy source terms models

SourceFiles
    macroModel.C
    newmacroModel.C

\*---------------------------------------------------------------------------*/

#ifndef macroModel_H
#define macroModel_H

#include "typeInfo.H"
#include "runTimeSelectionTables.H"
#include "volFields.H"
#include "dimensionedScalar.H"
#include "autoPtr.H"
#include "Pair.H"

#include "phaseModel.H"
#include "fvMatricesFwd.H"
#include "fvc.H"

#include "singleComponentSatProp.H"
#include "compressibleInterPhaseTransportModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class macroModel Declaration
\*---------------------------------------------------------------------------*/

class macroModel
    : public dictionary
{
private:
    //-
        dictionary macroModelCoeffs_;

protected:

        //- Thermo-package of phase 1
        const phaseModel& phase1_;

        //- Thermo-package of phase 2
        const phaseModel& phase2_;

        //- reference to pressure
        const volScalarField& p_;

        //- reference to saturation model
        singleComponentSatProp& satModel_;

        //- reference to turbulence Model
        const compressibleInterPhaseTransportModel& turbModel_;


    // Private Member Functions

        //- Disallow copy construct
        macroModel(const macroModel&);

        //- Disallow default bitwise assignment
        void operator=(const macroModel&);

public:

    //-
    dictionary& modelDict();

    //-
    const dictionary& modelDict() const;

public:

    //- Runtime type information
    TypeName("macroModel");

    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            macroModel,
            components,
            (
                const phaseModel& phase1,
                const phaseModel& phase2,
                const volScalarField& p,
                singleComponentSatProp& satModel,
                const compressibleInterPhaseTransportModel& turbModel,
                const dictionary& dict
            ),
            (phase1,phase2,p,satModel,turbModel,dict)
        );

    // Selectors

        //- Return a reference to the selected phaseChange model
        static autoPtr<macroModel> New
        (
            word macroModelTypeName,
            const phaseModel& phase1,
            const phaseModel& phase2,
            const volScalarField& p,
            singleComponentSatProp& satModel,
            const compressibleInterPhaseTransportModel& turbModel,
            const dictionary& dict
        );

    // Constructors

        //- Construct from components
        macroModel
        (
            const word& type,
            const phaseModel& phase1,
            const phaseModel& phase2,
            const volScalarField& p,
            singleComponentSatProp& satModel,
            const compressibleInterPhaseTransportModel& turbModel,
            const dictionary& dict
        );

    //- Destructor
    virtual ~macroModel()
    {}

    // Member Functions

        virtual void TSource1(fvScalarMatrix& T1Eqn) = 0;

        virtual void TSource2(fvScalarMatrix& T2Eqn) = 0;

        virtual void energySource(volScalarField& Q) = 0;

        virtual void energySource1(volScalarField& q1) = 0;

        virtual void energySource2(volScalarField& q2) = 0;

        virtual void massSource(volScalarField& rhoSource) = 0;

        virtual void alphaSource(volScalarField& rhoSource) = 0;


};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
