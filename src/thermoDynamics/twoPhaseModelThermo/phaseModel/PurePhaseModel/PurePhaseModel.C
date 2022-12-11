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

\*---------------------------------------------------------------------------*/

#include "PurePhaseModel.H"
#include "basicThermo.H"
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasePhaseModel, class phaseThermo>
Foam::PurePhaseModel<BasePhaseModel, phaseThermo>::PurePhaseModel
(
    const fvMesh& mesh,
    const dictionary& dict,
    const word& phaseName
)
:
    BasePhaseModel(mesh,dict, phaseName)
{
    thermoPtr_.reset
    (
        phaseThermo::New
        (
            mesh,
            phaseName
        ).ptr()
    );

}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasePhaseModel, class phaseThermo>
Foam::tmp<Foam::fvScalarMatrix>
Foam::PurePhaseModel<BasePhaseModel, phaseThermo>::YiEqn
(
    volScalarField& Yi,
    const volScalarField& muEff
)
{
    NotImplemented;
    return nullptr;
}


template<class BasePhaseModel, class phaseThermo>
const Foam::PtrList<Foam::volScalarField>&
Foam::PurePhaseModel<BasePhaseModel, phaseThermo>::Y() const
{
    return Y_;
}


template<class BasePhaseModel, class phaseThermo>
Foam::PtrList<Foam::volScalarField>&
Foam::PurePhaseModel<BasePhaseModel, phaseThermo>::Y()
{
    return Y_;
}


template<class BasePhaseModel, class phaseThermo>
const phaseThermo& Foam::PurePhaseModel<BasePhaseModel, phaseThermo>::
thermo() const
{
    return thermoPtr_();
}


template<class BasePhaseModel, class phaseThermo>
phaseThermo& Foam::PurePhaseModel<BasePhaseModel, phaseThermo>::thermo()
{
    return thermoPtr_();
}

template<class BasePhaseModel, class phaseThermo>
Foam::label Foam::PurePhaseModel<BasePhaseModel, phaseThermo>::inertIndex() const
{
    return -1;
}


// ************************************************************************* //
