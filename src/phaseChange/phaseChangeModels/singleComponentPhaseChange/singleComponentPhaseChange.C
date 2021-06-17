/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013 OpenFOAM Foundation
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


#include "fvCFD.H"
#include "singleComponentPhaseChange.H"


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::singleComponentPhaseChange::singleComponentPhaseChange
(
    const twoPhaseModelThermo& mixture,
    const volScalarField& p,
    const compressibleInterPhaseTransportModel& turbModel,
    reconstructionSchemes& surf
)
:
    IOdictionary
    (
        IOobject
        (
            "singleComponentPhaseChange",
            p.time().constant(),
            p.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        )
    ),
    // General data
    mesh_(mixture.phase1().mesh()),
    energyModel_(nullptr),
    macroModels_(),
    massModel_(nullptr),
    satProp_(nullptr),
    phase1_(mixture.phase1()),
    phase2_(mixture.phase2()),
    p_(p),
    turbModel_(turbModel),
    surf_(surf),
    psi0_
    (
        IOobject
        (
            "psi0_",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("0", dimDensity/dimTime, 0),
        "zeroGradient"
    ),
    massSource_
    (
        IOobject
        (
        "massSource_",
        mesh_.time().timeName(),
        mesh_,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
        mesh_,
        dimensionedScalar("0", dimDensity/dimTime, 0),
        "zeroGradient"
    ),
    alphaSource_
    (
        IOobject
        (
        "alphaSource_",
        mesh_.time().timeName(),
        mesh_,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
        mesh_,
        dimensionedScalar("0", dimless/dimTime, 0),
        "zeroGradient"
    ),
    limitHeatFlux_(false)
{
    IOdictionary phaseChangeProperties
    (
        IOobject
        (
            "phaseChangeProperties",
            mesh_.time().constant(),
            mesh_,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    );

    dictionary satPropertiesDict =
        phaseChangeProperties.subDict("satProperties");
    Info << "creating models" << endl;

    Info << "creating singleComponentSatProp Model" << endl;
    satProp_ =  singleComponentSatProp::New(mesh_,satPropertiesDict);

    Info << "creating Energy Source Term Model" << endl;
    energyModel_ =  energySourceTermModel::New
    (
        phase1_,
        phase2_,
        turbModel_,
        p_,
        satProp_.ref(),
        surf_,
        phaseChangeProperties
    );

    const dictionary marcoModelDict =
        phaseChangeProperties.subOrEmptyDict("marcoModels");

    label modeli = 0;
    macroModels_.setSize
    (
        marcoModelDict.size()
    );

    forAllConstIters(marcoModelDict, iter)
    {
        const word& key = iter().keyword();

        if (!marcoModelDict.isDict(key))
        {
            FatalErrorInFunction
                << "Found non-dictionary entry " << iter()
                << " in top-level dictionary " << marcoModelDict
                << exit(FatalError);
        }

        const dictionary& compmarcoModel = marcoModelDict.subDict(key);

        macroModels_.set
        (
            modeli,
            macroModel::New
            (
                word(compmarcoModel.lookup("type")),
                phase1_,
                phase2_,
                p_,
                satProp_.ref(),
                turbModel_,
                compmarcoModel
            )
        );

        ++modeli;
    }

    Info << "creating Mass Source Term Model" << endl;
    massModel_ =  massSourceTermModel::New
    (
        phase1_,
        phase2_,
        p_,
        satProp_.ref(),
        surf_,
        phaseChangeProperties
    );

}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::singleComponentPhaseChange::~singleComponentPhaseChange()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::singleComponentPhaseChange::correct()
{
    surf_.reconstruct(false);

    tmp<volScalarField> tphaseChangeEnergy = energyModel_->energySource();
    volScalarField& phaseChangeEnergy = tphaseChangeEnergy.ref();

    for (auto& mModel: macroModels_)
    {
        mModel.energySource(phaseChangeEnergy);
    }

    psi0_.ref() = phaseChangeEnergy.internalField() /
                  satProp_->L().internalField();



    const volScalarField& rho1 = phase1_.thermo().rho();
    if(limitHeatFlux_)
    {
        forAll(psi0_,celli)
        {
            scalar rhobyDt = rho1[celli]/mesh_.time().deltaTValue();
            scalar maxEvap = phase1_[celli]*rhobyDt; // positive
            scalar maxCond = -phase2_[celli]*rhobyDt; // negative
            psi0_[celli] = min(max(psi0_[celli],maxCond),maxEvap);
        }
    }
    psi0_.correctBoundaryConditions();

    massSource_ = massModel_->massSource(psi0_);

    for (auto& mModel: macroModels_)
    {
        mModel.massSource(massSource_);
    }

    alphaSource_ = massModel_->alphaSource(psi0_);

    for (auto& mModel: macroModels_)
    {
        mModel.alphaSource(alphaSource_);
    }


}

void Foam::singleComponentPhaseChange::correctSatProperties
(
    const volScalarField& p,
    const volScalarField& T
)
{
   satProp_->correct(p,T);
}


// ************************************************************************* //
