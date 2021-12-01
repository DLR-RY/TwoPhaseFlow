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

\*---------------------------------------------------------------------------*/


#include "Schrage.H"
#include "addToRunTimeSelectionTable.H"
#include "zeroGradientFvPatchFields.H"

#include "mathematicalConstants.H"
#include "surfaceInterpolate.H"
#include "fvcDiv.H"
#include "fvcGrad.H"



// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(Schrage, 0);
    addToRunTimeSelectionTable(energySourceTermModel,Schrage, components);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::Schrage::Schrage
(
    const phaseModel& phase1,
    const phaseModel& phase2,
    const compressibleInterPhaseTransportModel& turbModel,
    const volScalarField& p,
    singleComponentSatProp& satModel,
    reconstructionSchemes& surf,
    const dictionary& dict
)
:
    energySourceTermModel
    (
        typeName,
        phase1,
        phase2,
        turbModel,
        p,
        satModel,
        surf,
        dict
    ),
    evapCoeff_(modelDict().lookupOrDefault<scalar>("sigma",1)),
    R_("R",dimGasConstant,0)
{
    if (phase2_.thermo().incompressible())
    {
        R_.value() = modelDict().get<scalar>("R");
    }
}

// * * * * * * * * * * * * * * Protected Access Member Functions  * * * * *  //

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //


// ************************************************************************* //

Foam::tmp<Foam::fvScalarMatrix> Foam::Schrage::TSource1()
{
    surf_.reconstruct(false);
    const fvMesh& mesh = phase1_.mesh();

    const volScalarField& TSat = satModel_.TSat();

    volScalarField psi2 = phase2_.thermo().psi();
    if (phase2_.thermo().incompressible())
    {
        psi2 = 1/(R_*phase2_.thermo().T());
    }
    const volScalarField& rho2 = phase2_.thermo().rho();
    const volScalarField& T2 = phase2_.thermo().T();
    const volScalarField& T1 = phase1_.thermo().T();


    volScalarField TSource
    (
        (2*evapCoeff_)/(2-evapCoeff_)
        *pow(satModel_.L(),2)/(pow(2*constant::mathematical::pi/(psi2*T2),0.5))
        *rho2/pow(TSat,1.5)
    );

    TSource.ref() *= mag(surf_.normal().internalField())/mesh.V();

    tmp<fvScalarMatrix> T1Source(fvm::Sp(TSource,T1) - TSource*TSat);

    return T1Source;
}

Foam::tmp<Foam::fvScalarMatrix> Foam::Schrage::TSource2()
{
    surf_.reconstruct(false);
    const fvMesh& mesh = phase1_.mesh();

    const volScalarField& TSat = satModel_.TSat();

    volScalarField psi2 = phase2_.thermo().psi();
    if (phase2_.thermo().incompressible())
    {
        psi2 = 1/(R_*phase2_.thermo().T());
    }
    const volScalarField& rho2 = phase2_.thermo().rho();
    const volScalarField& T2 = phase2_.thermo().T();

    volScalarField TSource
    (
        (2*evapCoeff_)/(2-evapCoeff_)
        *pow(satModel_.L(),2)/(pow(2*constant::mathematical::pi/(psi2*T2),0.5))
        *rho2/pow(TSat,1.5)
    );

    TSource.ref() *= mag(surf_.normal().internalField())/mesh.V();

    tmp<fvScalarMatrix> T2Source(fvm::Sp(TSource,T2) - TSource*TSat);

    return T2Source;
}

Foam::tmp<Foam::volScalarField> Foam::Schrage::energySource()
{
    surf_.reconstruct(false);
    const volScalarField& TSat = satModel_.TSat();

    volScalarField psi2 = phase2_.thermo().psi();
    if (phase2_.thermo().incompressible())
    {
        psi2 = 1/(R_*phase2_.thermo().T());
    }
    const volScalarField& rho2 = phase2_.thermo().rho();
    const volScalarField& T1 = phase1_.thermo().T();
    const volScalarField& T2 = phase2_.thermo().T();

    volScalarField TSource
    (
        (2*evapCoeff_)/(2-evapCoeff_)
        *pow(satModel_.L(),2)/(pow(2*constant::mathematical::pi/(psi2*T2),0.5))
        *rho2/pow(TSat,1.5)
    );

    tmp<volScalarField> energySource(TSource*(T1-TSat) + TSource*(T2-TSat));
    volScalarField& energySourceRef = energySource.ref();
    energySourceRef.ref() *= mag(surf_.normal().internalField())/TSat.mesh().V();
    energySourceRef.correctBoundaryConditions();

    return energySource;

}


Foam::tmp<Foam::volScalarField> Foam::Schrage::energyFlux1()
{
    const volScalarField& TSat = satModel_.TSat();

    volScalarField psi2 = phase2_.thermo().psi();
    if (phase2_.thermo().incompressible())
    {
        psi2 = 1/(R_*phase2_.thermo().T());
    }
    const volScalarField& rho2 = phase2_.thermo().rho();
    const volScalarField& T1 = phase1_.thermo().T();
    const volScalarField& T2 = phase2_.thermo().T();

    volScalarField TSource
    (
        (2*evapCoeff_)/(2-evapCoeff_)
        *pow(satModel_.L(),2)/(pow(2*constant::mathematical::pi/(psi2*T2),0.5))
        *rho2/pow(TSat,1.5)
    );

    volScalarField interface(phase2_*0);
    interface.boundaryFieldRef() = Zero;

    for (const label celli: surf_.interfaceLabels())
    {
        interface[celli] = 1;
    }

    tmp<volScalarField> energyFlux1(TSource*(T1-TSat)*interface);

    return energyFlux1;

}

Foam::tmp<Foam::volScalarField> Foam::Schrage::energyFlux2()
{
    const volScalarField& TSat = satModel_.TSat();

    volScalarField psi2 = phase2_.thermo().psi();
    if (phase2_.thermo().incompressible())
    {
        psi2 = 1/(R_*phase2_.thermo().T());
    }
    const volScalarField& rho2 = phase2_.thermo().rho();
    const volScalarField& T2 = phase2_.thermo().T();

    volScalarField TSource
    (
        (2*evapCoeff_)/(2-evapCoeff_)
        *pow(satModel_.L(),2)/(pow(2*constant::mathematical::pi/(psi2*T2),0.5))
        *rho2/pow(TSat,1.5)
    );

    volScalarField interface(phase2_*0);
    interface.boundaryFieldRef() = Zero;

    for (const label celli: surf_.interfaceLabels())
    {
        interface[celli] = 1;
    }

    tmp<volScalarField> energyFlux2(TSource*(T2-TSat)*interface);

    return energyFlux2;

}
