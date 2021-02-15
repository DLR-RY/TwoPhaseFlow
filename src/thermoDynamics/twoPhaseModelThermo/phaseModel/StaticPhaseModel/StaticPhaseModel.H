/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2017 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

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
    Foam::StaticPhaseModel

Description
    Class which represents a static fluid phase.

SourceFiles
    StaticPhaseModel.C

\*---------------------------------------------------------------------------*/

#ifndef StaticPhaseModel_H
#define StaticPhaseModel_H

#include "phaseModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class phaseModel Declaration
\*---------------------------------------------------------------------------*/

template<class BasePhaseModel>
class StaticPhaseModel
:
    public BasePhaseModel
{
    // Private data

        //- Reference to U
        const volVectorField& U_;

        //- Flux
        surfaceScalarField phi_;

        //- Volumetric flux
        surfaceScalarField alphaPhi_;


public:

    // Constructors

        StaticPhaseModel
        (
            const fvMesh& mesh,
            const dictionary& dict,
            const word& phaseName
        );


    //- Destructor
    virtual ~StaticPhaseModel() = default;


    // Member Functions

        //- Correct the phase properties other than the thermo and turbulence
        virtual void correct();


        // Momentum

            //- Constant access the volumetric flux. Return zero field
            virtual tmp<surfaceScalarField> phi() const;

            //- Access the volumetric flux
            virtual const surfaceScalarField& phi();

            //- Constant access the volumetric flux of the phase.
            virtual tmp<surfaceScalarField> alphaPhi() const;

            //- Access the volumetric flux of the phase
            virtual surfaceScalarField& alphaPhi();

            //- Access const reference to U
            virtual tmp<volVectorField> U() const;

            //- Maximum diffusion number
            virtual tmp<surfaceScalarField> diffNo() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "StaticPhaseModel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
