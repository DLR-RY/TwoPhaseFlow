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


#include "constantHeatFlux.H"
#include "addToRunTimeSelectionTable.H"
#include "zeroGradientFvPatchFields.H"

#include "alphaContactAngleFvPatchScalarField.H"
#include "mathematicalConstants.H"
#include "surfaceInterpolate.H"
#include "fvcDiv.H"
#include "fvcGrad.H"



// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(constantHeatFlux, 0);
    addToRunTimeSelectionTable(energySourceTermModel,constantHeatFlux, components);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::constantHeatFlux::constantHeatFlux
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
    h_
    (
        IOobject
        (
            "h",
            phase2.mesh().time().timeName(),
            phase2.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        phase2.mesh(),
        dimensionedScalar
        (
            "h",
            dimPower/dimArea/dimTemperature,
            modelDict().get<scalar>("h")
        )
    )
{

}

// * * * * * * * * * * * * * * Protected Access Member Functions  * * * * *  //

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //


// ************************************************************************* //

Foam::tmp<Foam::fvScalarMatrix> Foam::constantHeatFlux::TSource1()
{
    surf_.reconstruct(false);
    const fvMesh& mesh = phase1_.mesh();
    
    const volScalarField& TSat = satModel_.TSat();

    const volScalarField& T1 = phase1_.thermo().T();

    volScalarField TSource = h_;

    TSource.ref() *= mag(surf_.normal().internalField())/mesh.V();

    tmp<fvScalarMatrix> T1Source(fvm::Sp(TSource,T1) - TSource*TSat);

    return T1Source;
}

Foam::tmp<Foam::fvScalarMatrix> Foam::constantHeatFlux::TSource2()
{
    surf_.reconstruct(false);
    const fvMesh& mesh = phase1_.mesh();

    const volScalarField& TSat = satModel_.TSat();

    const volScalarField& T2 = phase2_.thermo().T();

    volScalarField TSource = h_;

    TSource.ref() *= mag(surf_.normal().internalField())/mesh.V();

    tmp<fvScalarMatrix> T2Source(fvm::Sp(TSource,T2) - TSource*TSat);

    return T2Source;
}

Foam::tmp<Foam::volScalarField> Foam::constantHeatFlux::energySource()
{
    surf_.reconstruct(false);
    const volScalarField& TSat = satModel_.TSat();

    const volScalarField& T1 = phase1_.thermo().T();
    const volScalarField& T2 = phase2_.thermo().T();

    volScalarField TSource = h_;

    tmp<volScalarField> energySource(TSource*(T1-TSat) + TSource*(T2-TSat));
    volScalarField& energySourceRef = energySource.ref();
    energySourceRef.ref() *= mag(surf_.normal().internalField())/TSat.mesh().V();
    energySourceRef.correctBoundaryConditions();

    return energySource;

}


Foam::tmp<Foam::volScalarField> Foam::constantHeatFlux::energyFlux1()
{
    const volScalarField& TSat = satModel_.TSat();

    const volScalarField& T1 = phase1_.thermo().T();

    volScalarField TSource = h_;

    volScalarField interface = phase2_*0;
    interface.boundaryFieldRef() = Zero;

    for (const label celli: surf_.interfaceLabels())
    {
        interface[celli] = 1;
    }

    tmp<volScalarField> energyFlux1(TSource*(T1-TSat)*interface);

    return energyFlux1;

}

Foam::tmp<Foam::volScalarField> Foam::constantHeatFlux::energyFlux2()
{
    const volScalarField& TSat = satModel_.TSat();

    const volScalarField& T2 = phase2_.thermo().T();

    volScalarField TSource = h_;

    volScalarField interface = phase2_*0;
    interface.boundaryFieldRef() = Zero;

    for (const label celli: surf_.interfaceLabels())
    {
        interface[celli] = 1;
    }

    tmp<volScalarField> energyFlux2(TSource*(T2-TSat)*interface);

    return energyFlux2;

}
