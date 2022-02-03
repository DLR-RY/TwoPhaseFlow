/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013-2017 OpenFOAM Foundation
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

#include "twoPhaseModelThermo.H"
#include "gradientEnergyFvPatchScalarField.H"
#include "mixedEnergyFvPatchScalarField.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(twoPhaseModelThermo, 0);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::twoPhaseModelThermo::twoPhaseModelThermo
(
    const volVectorField& U,
    const surfaceScalarField& phi
)
:
    psiThermo(U.mesh(), word::null),
    phase1Name_(wordList(lookup("phases"))[0]),
    phase2Name_(wordList(lookup("phases"))[1]),
    phase1_(phaseModel::New(U.mesh(), *this, phase1Name_)),
    phase2_(phaseModel::New(U.mesh(), *this, phase2Name_)),
    he_
    (
        IOobject
        (
            "twoPhaseModelThermo::he",
            U.mesh().time().timeName(),
            U.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        U.mesh(),
        dimensionedScalar("0", dimEnergy/dimMass,0)
    )
{
    alpha2() = scalar(1) - alpha1();

    phase1_->thermo().validate(phase1Name_, "e");
    phase2_->thermo().validate(phase2Name_, "e");

    correctThermo(phase1_->thermo().T(),phase2_->thermo().T());
    correct();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::twoPhaseModelThermo::~twoPhaseModelThermo()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::twoPhaseModelThermo::nearInterface() const
{
    return pos0(phase1_() - 0.01)*pos0(0.99 - phase1_());
}


void Foam::twoPhaseModelThermo::correctThermo()
{
    phase1_->thermo().he() = phase1_->thermo().he(p_, T_);
    phase1_->correct();

    phase2_->thermo().he() = phase2_->thermo().he(p_, T_);
    phase2_->correct();

    he_ = phase1_()*phase1_->thermo().he() + phase2_()*phase2_->thermo().he();
}


void Foam::twoPhaseModelThermo::correctThermo
(
      const volScalarField& TL,
      const volScalarField& TV
)
{
    phase1_->thermo().he() = phase1_->thermo().he(p_, TL);
    phase1_->correct();

    phase2_->thermo().he() = phase2_->thermo().he(p_, TV);
    phase2_->correct();
    // wallHeatFlux uses alpha = kappa/Cp instead of kappa/Cv therefore
    // he needs the be scaled
    he_ = phase1_()*phase1_->thermo().he()*phase1_->thermo().CpByCpv()
        + phase2_()*phase2_->thermo().he()*phase2_->thermo().CpByCpv();
}


void Foam::twoPhaseModelThermo::correct()
{
    psi_ =
        phase1_()*phase1_->thermo().psi() + phase2_()*phase2_->thermo().psi();
    mu_ =
        phase1_()*phase1_->thermo().mu() + phase2_()*phase2_->thermo().mu();
    alpha_ =
        phase1_()*phase1_->thermo().alpha()
      + phase2_()*phase2_->thermo().alpha();

}

Foam::word Foam::twoPhaseModelThermo::thermoName() const
{
    return
        phase1_->thermo().thermoName() + ',' + phase2_->thermo().thermoName();
}

bool Foam::twoPhaseModelThermo::incompressible() const
{
    return phase1_->thermo().incompressible()
        && phase2_->thermo().incompressible();
}


bool Foam::twoPhaseModelThermo::isochoric() const
{
    return phase1_->thermo().isochoric() && phase2_->thermo().isochoric();
}


Foam::tmp<Foam::volScalarField> Foam::twoPhaseModelThermo::he
(
    const volScalarField& p,
    const volScalarField& T
) const
{
    return phase1_()*phase1_->thermo().he(p, T)
         + phase2_()*phase2_->thermo().he(p, T);
}


Foam::tmp<Foam::scalarField> Foam::twoPhaseModelThermo::he
(
    const scalarField& p,
    const scalarField& T,
    const labelList& cells
) const
{
    return
        scalarField(phase1_, cells)*phase1_->thermo().he(p, T, cells)
      + scalarField(phase2_, cells)*phase2_->thermo().he(p, T, cells);
}


Foam::tmp<Foam::scalarField> Foam::twoPhaseModelThermo::he
(
    const scalarField& p,
    const scalarField& T,
    const label patchi
) const
{
    return
        phase1_().boundaryField()[patchi]*phase1_->thermo().he(p, T, patchi)
      + phase2_().boundaryField()[patchi]*phase2_->thermo().he(p, T, patchi);
}


Foam::tmp<Foam::volScalarField> Foam::twoPhaseModelThermo::hc() const
{
    return phase1_()*phase1_->thermo().hc() + phase2_()*phase2_->thermo().hc();
}


Foam::tmp<Foam::scalarField> Foam::twoPhaseModelThermo::THE
(
    const scalarField& h,
    const scalarField& p,
    const scalarField& T0,
    const labelList& cells
) const
{
    NotImplemented;
    return T0;
}


Foam::tmp<Foam::scalarField> Foam::twoPhaseModelThermo::THE
(
    const scalarField& h,
    const scalarField& p,
    const scalarField& T0,
    const label patchi
) const
{
    NotImplemented;
    return T0;
}


Foam::tmp<Foam::volScalarField> Foam::twoPhaseModelThermo::Cp() const
{
    return phase1_()*phase1_->thermo().Cp() + phase2_()*phase2_->thermo().Cp();
}


Foam::tmp<Foam::scalarField> Foam::twoPhaseModelThermo::Cp
(
    const scalarField& p,
    const scalarField& T,
    const label patchi
) const
{
    return
        phase1_().boundaryField()[patchi]*phase1_->thermo().Cp(p, T, patchi)
      + phase2_().boundaryField()[patchi]*phase2_->thermo().Cp(p, T, patchi);
}


Foam::tmp<Foam::volScalarField> Foam::twoPhaseModelThermo::Cv() const
{
    return phase1_()*phase1_->thermo().Cv() + phase2_()*phase2_->thermo().Cv();
}


Foam::tmp<Foam::scalarField> Foam::twoPhaseModelThermo::Cv
(
    const scalarField& p,
    const scalarField& T,
    const label patchi
) const
{
    return
        phase1_().boundaryField()[patchi]*phase1_->thermo().Cv(p, T, patchi)
      + phase2_().boundaryField()[patchi]*phase2_->thermo().Cv(p, T, patchi);
}

Foam::tmp<Foam::scalarField> Foam::twoPhaseModelThermo::rhoEoS
(
    const scalarField& p,
    const scalarField& T,
    const labelList& cells
) const
{
    NotImplemented;
    return nullptr;
}


Foam::tmp<Foam::volScalarField> Foam::twoPhaseModelThermo::gamma() const
{
    return phase1_()*phase1_->thermo().gamma()
         + phase2_()*phase2_->thermo().gamma();
}


Foam::tmp<Foam::scalarField> Foam::twoPhaseModelThermo::gamma
(
    const scalarField& p,
    const scalarField& T,
    const label patchi
) const
{
    return
        phase1_().boundaryField()[patchi]*phase1_->thermo().gamma(p, T, patchi)
      + phase2_().boundaryField()[patchi]*phase2_->thermo().gamma(p, T, patchi);
}


Foam::tmp<Foam::volScalarField> Foam::twoPhaseModelThermo::Cpv() const
{
    return phase1_()*phase1_->thermo().Cpv()
         + phase2_()*phase2_->thermo().Cpv();
}


Foam::tmp<Foam::scalarField> Foam::twoPhaseModelThermo::Cpv
(
    const scalarField& p,
    const scalarField& T,
    const label patchi
) const
{
    return
        phase1_().boundaryField()[patchi]*phase1_->thermo().Cpv(p, T, patchi)
      + phase2_().boundaryField()[patchi]*phase2_->thermo().Cpv(p, T, patchi);
}


Foam::tmp<Foam::volScalarField> Foam::twoPhaseModelThermo::CpByCpv() const
{
    return
        phase1_()*phase1_->thermo().CpByCpv()
      + phase2_()*phase2_->thermo().CpByCpv();
}

Foam::tmp<Foam::volScalarField> Foam::twoPhaseModelThermo::W() const
{
    return phase1_()*phase1_->thermo().W() + phase2_()*phase1_->thermo().W();
}


Foam::tmp<Foam::scalarField> Foam::twoPhaseModelThermo::CpByCpv
(
    const scalarField& p,
    const scalarField& T,
    const label patchi
) const
{
    return
        phase1_().boundaryField()[patchi]*phase1_->thermo().CpByCpv(p, T, patchi)
      + phase2_().boundaryField()[patchi]*phase2_->thermo().CpByCpv(p, T, patchi);
}


Foam::tmp<Foam::volScalarField> Foam::twoPhaseModelThermo::nu() const
{
    return mu()/(phase1_()*phase1_->thermo().rho()
         + phase2_()*phase2_->thermo().rho());
}


Foam::tmp<Foam::scalarField> Foam::twoPhaseModelThermo::nu
(
    const label patchi
) const
{
    return
        mu(patchi)
       /(
            phase1_().boundaryField()[patchi]*phase1_->thermo().rho(patchi)
          + phase2_().boundaryField()[patchi]*phase2_->thermo().rho(patchi)
        );
}


Foam::tmp<Foam::volScalarField> Foam::twoPhaseModelThermo::kappa() const
{
    return phase1_()*phase1_->thermo().kappa()
         + phase2_()*phase2_->thermo().kappa();
}


Foam::tmp<Foam::scalarField> Foam::twoPhaseModelThermo::kappa
(
    const label patchi
) const
{
    return
        phase1_().boundaryField()[patchi]*phase1_->thermo().kappa(patchi)
      + phase2_().boundaryField()[patchi]*phase2_->thermo().kappa(patchi);
}


Foam::tmp<Foam::volScalarField> Foam::twoPhaseModelThermo::alphahe() const
{
    return phase1_()*phase1_->thermo().alphahe()
         + phase2_()*phase2_->thermo().alphahe();
}


Foam::tmp<Foam::scalarField> Foam::twoPhaseModelThermo::alphahe
(
    const label patchi
) const
{
    return
        phase1_().boundaryField()[patchi]*phase1_->thermo().alphahe(patchi)
      + phase2_().boundaryField()[patchi]*phase2_->thermo().alphahe(patchi);
}



Foam::tmp<Foam::volScalarField> Foam::twoPhaseModelThermo::kappaEff
(
    const volScalarField& alphat
) const
{
    return
        phase1_()*phase1_->thermo().kappaEff(alphat)
      + phase2_()*phase2_->thermo().kappaEff(alphat);
}


Foam::tmp<Foam::scalarField> Foam::twoPhaseModelThermo::kappaEff
(
    const scalarField& alphat,
    const label patchi
) const
{
    return
        phase1_().boundaryField()[patchi]*phase1_->thermo().kappaEff(alphat, patchi)
      + phase2_().boundaryField()[patchi]*phase2_->thermo().kappaEff(alphat, patchi);
}


Foam::tmp<Foam::volScalarField> Foam::twoPhaseModelThermo::alphaEff
(
    const volScalarField& alphat
) const
{
    return
        phase1_()*phase1_->thermo().alphaEff(alphat)
      + phase2_()*phase2_->thermo().alphaEff(alphat);
}


Foam::tmp<Foam::scalarField> Foam::twoPhaseModelThermo::alphaEff
(
    const scalarField& alphat,
    const label patchi
) const
{
    return
        phase1_().boundaryField()[patchi]*phase1_->thermo().alphaEff(alphat, patchi)
      + phase2_().boundaryField()[patchi]*phase2_->thermo().alphaEff(alphat, patchi);
}


bool Foam::twoPhaseModelThermo::read()
{
    if (psiThermo::read())
    {
        return true; //interfaceProperties::read();
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
