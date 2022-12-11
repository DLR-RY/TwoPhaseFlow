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

#include "MultiComponentPhaseModel.H"



#include "fvmDdt.H"
#include "fvmDiv.H"
#include "fvmSup.H"
#include "fvmLaplacian.H"
#include "fvcDdt.H"
#include "fvcDiv.H"
#include "fvcDDt.H"
#include "fvMatrix.H"
#include "fvcFlux.H"
#include "CMULES.H"
#include "subCycle.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasePhaseModel, class phaseThermo>
Foam::MultiComponentPhaseModel<BasePhaseModel, phaseThermo>::
MultiComponentPhaseModel
(
    const fvMesh& mesh,
    const dictionary& dict,
    const word& phaseName
)
:
    BasePhaseModel(mesh, dict, phaseName),
    species_(),
    inertIndex_(-1),
    Sc_
    (
        "Sc",
        dimless,
        dict.subDict(phaseName)
    ),
    residualAlpha_
    (
        "residualAlpha",
        dimless,
        mesh.solverDict("Yi")
    )
{
    thermoPtr_.reset
    (
        phaseThermo::New
        (
            mesh,
            phaseName
        ).ptr()
    );

    const word inertSpecie
    (
        thermoPtr_->lookupOrDefault("inertSpecie", word::null)
    );

    if (inertSpecie != word::null)
    {
        inertIndex_ = thermoPtr_->composition().species()[inertSpecie];
    }

    if (thermoPtr_->composition().species().size() == 0)
    {
        FatalErrorInFunction
            << " The selected thermo is pure. Use a multicomponent thermo."
            << exit(FatalError);
    }

    species_ = thermoPtr_->composition().species();
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasePhaseModel, class phaseThermo>
const phaseThermo&
Foam::MultiComponentPhaseModel<BasePhaseModel, phaseThermo>::thermo() const
{
    return thermoPtr_();
}


template<class BasePhaseModel, class phaseThermo>
phaseThermo&
Foam::MultiComponentPhaseModel<BasePhaseModel, phaseThermo>::thermo()
{
    return thermoPtr_();
}

template<class BasePhaseModel, class phaseThermo>
void Foam::MultiComponentPhaseModel<BasePhaseModel, phaseThermo>::correct()
{
    volScalarField Yt
    (
        IOobject
        (
            IOobject::groupName("Yt", this->name()),
            this->mesh().time().timeName(),
            this->mesh()
        ),
        this->mesh(),
        dimensionedScalar(dimless, Zero)
    );

    PtrList<volScalarField>& Yi = Y();

    forAll(Yi, i)
    {
        if (i != inertIndex_)
        {
            Yt += Yi[i];
        }
    }

    if (inertIndex_ != -1)
    {
        Yi[inertIndex_] = scalar(1) - Yt;
        Yi[inertIndex_].max(0);
    }
    else
    {
        forAll(Yi, i)
        {
            Yi[i] /= Yt;
            Yi[i].max(0);
        }
    }

    return thermo().correct();
}


template<class BasePhaseModel, class phaseThermo>
Foam::tmp<Foam::fvScalarMatrix>
Foam::MultiComponentPhaseModel<BasePhaseModel, phaseThermo>::YiEqn
(
    volScalarField& Yi,
    const volScalarField& muEff
)
{
    if
    (
        (inertIndex_ != -1)
     && (
            (
                Yi.name()
             == IOobject::groupName
                (
                    this->thermoPtr_->composition().species()[inertIndex_],
                    this->name()
                )
            )
         || (
               !this->thermoPtr_->composition().active
                (
                    this->thermoPtr_->composition().species()[Yi.member()]
                )
            )
        )
    )
    {
        return nullptr;
    }

    const volScalarField& alpha = *this;
    const surfaceScalarField alphaRhoPhi(this->alphaPhi()*fvc::interpolate(this->rho()));
    const volScalarField& rho = this->rho();
    const volScalarField continuityError(fvc::ddt(alpha, rho) + fvc::div(alphaRhoPhi));

    return
    (
        fvm::ddt(alpha, rho, Yi)
      + fvm::div(alphaRhoPhi, Yi, "div(" + alphaRhoPhi.name() + ",Yi)")
      - fvm::Sp(continuityError, Yi)

      - fvm::laplacian
        (
            fvc::interpolate(alpha)
           *fvc::interpolate(muEff/Sc_),
            Yi
        )
     ==
        fvc::ddt(residualAlpha_*rho, Yi)
      - fvm::ddt(residualAlpha_*rho, Yi)
    );
}


template<class BasePhaseModel, class phaseThermo>
const Foam::PtrList<Foam::volScalarField>&
Foam::MultiComponentPhaseModel<BasePhaseModel, phaseThermo>::Y() const
{
    return thermoPtr_->composition().Y();
}


template<class BasePhaseModel, class phaseThermo>
Foam::PtrList<Foam::volScalarField>&
Foam::MultiComponentPhaseModel<BasePhaseModel, phaseThermo>::Y()
{
    return thermoPtr_->composition().Y();
}


template<class BasePhaseModel, class phaseThermo>
Foam::label
Foam::MultiComponentPhaseModel<BasePhaseModel, phaseThermo>::inertIndex() const
{
    return inertIndex_;
}


// ************************************************************************* //
