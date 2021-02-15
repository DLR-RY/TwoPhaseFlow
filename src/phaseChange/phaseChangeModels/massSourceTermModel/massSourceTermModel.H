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
    Foam::massSourceTermModel

Description
    base class of the mass source term model

SourceFiles
    massSourceTermModel.C
    newMassSourceTermModel.C

\*---------------------------------------------------------------------------*/

#ifndef massSourceTermModel_H
#define massSourceTermModel_H

#include "typeInfo.H"
#include "runTimeSelectionTables.H"
#include "volFields.H"
#include "dimensionedScalar.H"
#include "autoPtr.H"
#include "Pair.H"

#include "phaseModel.H"
#include "fvMatricesFwd.H"
#include "singleComponentSatProp.H"
#include "reconstructionSchemes.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class massSourceTermModel Declaration
\*---------------------------------------------------------------------------*/

class massSourceTermModel
    : public dictionary
{

private:

    //-
        dictionary massSourceTermModelCoeffs_;

protected:

        //- Thermo-package of phase 1
        const phaseModel& phase1_;

        //- Thermo-package of phase 2
        const phaseModel& phase2_;

        const volScalarField& p_;

        singleComponentSatProp& satModel_;

        reconstructionSchemes& surf_;

    // Private Member Functions

        //- Disallow copy construct
        massSourceTermModel(const massSourceTermModel&);

        //- Disallow default bitwise assignment
        void operator=(const massSourceTermModel&);

public:

    //-
    dictionary& modelDict();

    //-
    const dictionary& modelDict() const;

public:

    //- Runtime type information
    TypeName("massSourceTermModel");

    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            massSourceTermModel,
            components,
            (
                const phaseModel& phase1,
                const phaseModel& phase2,
                const volScalarField& p,
                singleComponentSatProp& satModel,
                reconstructionSchemes& surf,
                const dictionary& dict
            ),
            (phase1,phase2,p,satModel,surf,dict)
        );

    // Selectors

        //- Return a reference to the selected phaseChange model
        static autoPtr<massSourceTermModel> New
        (
            const phaseModel& phase1,
            const phaseModel& phase2,
            const volScalarField& p,
            singleComponentSatProp& satModel,
            reconstructionSchemes& surf,
            const dictionary& dict
        );

    // Constructors

        //- Construct from components
        massSourceTermModel
        (
            const word& type,
            const phaseModel& phase1,
            const phaseModel& phase2,
            const volScalarField& p,
            singleComponentSatProp& satModel,
            reconstructionSchemes& surf,
            const dictionary& dict
        );

    //- Destructor
    virtual ~massSourceTermModel()
    {}

    // Member Functions

        virtual tmp<volScalarField> massSource(volScalarField& rhoSource) = 0;

        virtual tmp<volScalarField> alphaSource(volScalarField& rhoSource) = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
