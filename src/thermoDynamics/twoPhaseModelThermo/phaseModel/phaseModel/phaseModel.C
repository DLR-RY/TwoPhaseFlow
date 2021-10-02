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

#include "phaseModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(phaseModel, 0);
    defineRunTimeSelectionTable(phaseModel, dictionary);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::phaseModel::phaseModel
(
    const fvMesh& mesh,
    const dictionary& dict,
    const word& phaseName
)
:
    volScalarField
    (
        IOobject
        (
            IOobject::groupName("alpha", phaseName),
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("alpha", dimless, 0)
    ),
    mesh_(mesh),
    name_(phaseName)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::word& Foam::phaseModel::name() const
{
    return name_;
}



void Foam::phaseModel::correct()
{
    thermo().correct();
}


void Foam::phaseModel::correctTurbulence()
{
    // do nothing
}


Foam::tmp<Foam::volScalarField> Foam::phaseModel::rho() const
{
    return thermo().rho();
}


Foam::tmp<Foam::scalarField> Foam::phaseModel::rho(const label patchI) const
{
     return thermo().rho(patchI);
}


Foam::tmp<Foam::volScalarField> Foam::phaseModel::hc() const
{
     return thermo().hc();
}


Foam::tmp<Foam::volScalarField> Foam::phaseModel::Cp() const
{
     return thermo().Cp();
}


Foam::tmp<Foam::scalarField> Foam::phaseModel::Cp
(
    const scalarField& p,
    const scalarField& T,
    const label patchI
) const
{
    return (thermo().Cp(p, T, patchI));
}


Foam::tmp<Foam::volScalarField> Foam::phaseModel::Cv() const
{
    return thermo().Cv();
}


Foam::tmp<Foam::scalarField> Foam::phaseModel::Cv
(
    const scalarField& p,
    const scalarField& T,
    const label patchI
) const
{
    return thermo().Cv(p, T, patchI);
}


Foam::tmp<Foam::volScalarField> Foam::phaseModel::gamma() const
{
    return thermo().gamma();
}


Foam::tmp<Foam::scalarField> Foam::phaseModel::gamma
(
    const scalarField& p,
    const scalarField& T,
    const label patchI
) const
{
    return thermo().gamma(p, T, patchI);
}


Foam::tmp<Foam::volScalarField> Foam::phaseModel::Cpv() const
{
    return thermo().Cpv();
}


Foam::tmp<Foam::scalarField> Foam::phaseModel::Cpv
(
    const scalarField& p,
    const scalarField& T,
    const label patchI
) const
{
    return thermo().Cpv(p, T, patchI);
}


Foam::tmp<Foam::volScalarField> Foam::phaseModel::CpByCpv() const
{
     return thermo().CpByCpv();
}


Foam::tmp<Foam::scalarField> Foam::phaseModel::CpByCpv
(
    const scalarField& p,
    const scalarField& T,
    const label patchI
) const
{
    return thermo().CpByCpv(p, T, patchI);
}


const Foam::volScalarField& Foam::phaseModel::alpha() const
{
    return thermo().alpha();
}


const Foam::scalarField& Foam::phaseModel::alpha(const label patchI) const
{
    return thermo().alpha(patchI);
}


Foam::tmp<Foam::volScalarField> Foam::phaseModel::kappa() const
{
    return thermo().kappa();
}


Foam::tmp<Foam::scalarField> Foam::phaseModel::kappa(const label patchI) const
{
    return thermo().kappa(patchI);
}

Foam::tmp<Foam::volScalarField> Foam::phaseModel::alphahe() const
{
    return thermo().alphahe();
}


Foam::tmp<Foam::scalarField> Foam::phaseModel::alphahe(const label patchI) const
{
    return thermo().alphahe(patchI);
}


Foam::tmp<Foam::volScalarField>Foam::phaseModel::kappaEff
(
    const volScalarField& kappat
) const
{
    tmp<volScalarField> kappaEff(kappa() + kappat);
    kappaEff.ref().rename("kappaEff" + name_);
    return kappaEff;
}


Foam::tmp<Foam::scalarField> Foam::phaseModel::kappaEff
(
    const scalarField& kappat,
    const label patchI
) const
{
    return (kappa(patchI) + kappat);
}


Foam::tmp<Foam::volScalarField> Foam::phaseModel::alphaEff
(
    const volScalarField& alphat
) const
{
    return (thermo().alpha() + alphat);
}


Foam::tmp<Foam::scalarField> Foam::phaseModel::alphaEff
(
    const scalarField& alphat,
    const label patchI
) const
{
    return (thermo().alpha(patchI) + alphat);
}


Foam::tmp<Foam::volScalarField> Foam::phaseModel::mu() const
{
    return thermo().mu();
}


Foam::tmp<Foam::scalarField> Foam::phaseModel::mu(const label patchi) const
{
    return thermo().mu(patchi);
}


Foam::tmp<Foam::volScalarField> Foam::phaseModel::nu() const
{
    return thermo().nu();
}


Foam::tmp<Foam::scalarField> Foam::phaseModel::nu(const label patchi) const
{
    return thermo().nu(patchi);
}


bool Foam::phaseModel::read()
{
    return true;
}


// ************************************************************************* //
