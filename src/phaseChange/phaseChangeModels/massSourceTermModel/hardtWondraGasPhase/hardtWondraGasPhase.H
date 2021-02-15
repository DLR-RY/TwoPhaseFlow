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
    Foam::hardtWondraGasPhase

Description
    Smears the mass source terms according to:

    S. Hardt, F. Wondra,
    Evaporation model for interfacial flows based on a
    continuum-field representation of the source terms,
    Journal of Computational Physics 227 (11) (2008) 5871–5895
    doi:10.1016/j.jcp.2008.02.020

    but mass source in the gas phase are applied in the VoF equation

SourceFiles
    hardtWondraGasPhase.C

\*---------------------------------------------------------------------------*/

#ifndef hardtWondraGasPhase_H
#define hardtWondraGasPhase_H

#include "typeInfo.H"
#include "massSourceTermModel.H"
#include "volFields.H"
#include "dimensionedScalar.H"
#include "autoPtr.H"

#include "fvScalarMatrix.H"
#include "fvm.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class hardtWondraGasPhase Declaration
\*---------------------------------------------------------------------------*/

class hardtWondraGasPhase
    : public massSourceTermModel
{
private:

    scalar cutoff_;
    scalar spread_;

    // Private Member Functions

        //- Disallow copy construct
        hardtWondraGasPhase(const hardtWondraGasPhase&);

        //- Disallow default bitwise assignment
        void operator=(const hardtWondraGasPhase&);

public:

    //- Runtime type information
    TypeName("hardtWondraGasPhase");

    // Constructors

        //- Construct from components
        hardtWondraGasPhase
        (
            const phaseModel& phase1,
            const phaseModel& phase2,
            const volScalarField& p,
            singleComponentSatProp& satModel,
            reconstructionSchemes& surf,
            const dictionary& dict
        );

    //- Destructor
    virtual ~hardtWondraGasPhase()
    {}

    virtual tmp<volScalarField> massSource(volScalarField& rhoSource);

    virtual tmp<volScalarField> alphaSource(volScalarField& rhoSource);

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
