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


#include "directEvaporation.H"
#include "addToRunTimeSelectionTable.H"
#include "zeroGradientFvPatchFields.H"

#include "alphaContactAngleFvPatchScalarField.H"
#include "mathematicalConstants.H"
#include "surfaceInterpolate.H"
#include "fvcDiv.H"
#include "fvcGrad.H"
#include "fvcSnGrad.H"

#include "centredCPCCellToCellStencilObject.H"
#include "SortableList.H"

#include "plane.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(directEvaporation, 0);
    addToRunTimeSelectionTable(energySourceTermModel,directEvaporation, components);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::directEvaporation::directEvaporation
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
    impDiffFlux_(phase1.mesh()),
    superheated_(modelDict().get<scalar>("superheated"))
{
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //


// ************************************************************************* //

Foam::tmp<Foam::fvScalarMatrix> Foam::directEvaporation::TSource1()
{
    surf_.reconstruct(false);
    const volScalarField& TSat = satModel_.TSat();
    const volScalarField& T1 = phase1_.thermo().T();

    tmp<fvScalarMatrix> TSource1
    (
        impDiffFlux_.diffusiveFlux
        (
            surf_.interfaceCell(),
            surf_.normal(),
            surf_.centre(),
            phase1_.thermo().T(),
            phase1_.thermo().kappaEff(turbModel_.alphat()),
            TSat,
            false,
            dimPower
        )
    );
    fvScalarMatrix& TSource1Ref = TSource1.ref();

    const volVectorField& normal = surf_.normal();
    const fvMesh& mesh = T1.mesh();

    // assume that alpha1 is the more dense phase
    dimensionedScalar superHeated("superHeated",dimTemperature,superheated_);
    dimensionedScalar deltaT = mesh.time().deltaT();

    volScalarField mEvap_Boiling
    (
        "mEvap_Boiling",
        phase1_*pos0(T1.oldTime()-(TSat+superHeated))
        *(T1.oldTime()-(TSat+superHeated))
        *phase1_.thermo().Cp()*phase1_.thermo().rho()/deltaT
    );
    if (mesh.time().outputTime())
    {
        mEvap_Boiling.write();
    }

    mEvap_Boiling.correctBoundaryConditions();

    TSource1Ref.source() -= mEvap_Boiling.internalField()*mesh.V();

    return TSource1;
}


Foam::tmp<Foam::fvScalarMatrix> Foam::directEvaporation::TSource2()
{
    surf_.reconstruct(false);
    const volScalarField& TSat = satModel_.TSat();

    tmp<fvScalarMatrix> TSource2
    (
        impDiffFlux_.diffusiveFlux
        (
            surf_.interfaceCell(),
            surf_.normal(),
            surf_.centre(),
            phase2_.thermo().T(),
            phase2_.thermo().kappaEff(turbModel_.alphat()),
            TSat,
            true,
            dimPower
        )
    );

    return TSource2;
}


Foam::tmp<Foam::volScalarField> Foam::directEvaporation::energySource()
{
    surf_.reconstruct(false);
    tmp<volScalarField> energySource
    (
        energyFlux1() + energyFlux2()
    );

    volScalarField& energySourceRef = energySource.ref();
    energySourceRef.ref() *= mag(surf_.normal().internalField())/phase1_.mesh().V();
    energySourceRef.correctBoundaryConditions();

    return energySource;

}


Foam::tmp<Foam::volScalarField> Foam::directEvaporation::energyFlux1()
{
    const volScalarField& TSat = satModel_.TSat();
    const volScalarField& T1 = phase1_.thermo().T();
    surf_.reconstruct(false);

    tmp<volScalarField> energyFlux1
    (
        impDiffFlux_.diffusiveFlux
        (
            surf_.interfaceCell(),
            surf_.normal(),
            surf_.centre(),
            phase1_.thermo().T(),
            phase1_.thermo().kappaEff(turbModel_.alphat()),
            TSat,
            false
        )
    );
    volScalarField& energyFluxRef1 = energyFlux1.ref();

    const volVectorField& normal = surf_.normal();
    const fvMesh& mesh = T1.mesh();

    // assume that alpha1 is the more dense phase
    dimensionedScalar superHeated("superHeated",dimTemperature,superheated_);
    dimensionedScalar deltaT = mesh.time().deltaT();

    volScalarField mEvap_Boiling
    (
        "mEvap_Boiling",
        phase1_*pos0(T1.oldTime()-(TSat+superHeated))
        *(T1.oldTime()-(TSat+superHeated))
        *phase1_.thermo().Cp()*phase1_.thermo().rho()/deltaT
    );

    if (mesh.time().outputTime())
    {
        mEvap_Boiling.write();
    }

    mEvap_Boiling.correctBoundaryConditions();

    volScalarField interface
    (
        "interface",
        mag(normal)/(mag(normal)+dimensionedScalar("0",dimArea,1e-16))
    );
    dimensionedScalar evapMass = fvc::domainIntegrate(mEvap_Boiling);

    volScalarField scaledInterface = interface*dimensionedScalar("0",dimless,1);
    dimensionedScalar area("area",dimArea,gSum(mag(normal)().internalField()));
    if (area.value() != 0)
    {
        dimensionedScalar scale = evapMass/area;
        scaledInterface *= scale;
    }

    energyFluxRef1 += scaledInterface;

    return energyFlux1;
}


Foam::tmp<Foam::volScalarField> Foam::directEvaporation::energyFlux2()
{
    const volScalarField& TSat = satModel_.TSat();
    surf_.reconstruct(false);

    Foam::tmp<Foam::volScalarField> energyFlux2
    (
        impDiffFlux_.diffusiveFlux
        (
            surf_.interfaceCell(),
            surf_.normal(),
            surf_.centre(),
            phase2_.thermo().T(),
            phase2_.thermo().kappaEff(turbModel_.alphat()),
            TSat,
            true
        )
    );

    return energyFlux2;
}
