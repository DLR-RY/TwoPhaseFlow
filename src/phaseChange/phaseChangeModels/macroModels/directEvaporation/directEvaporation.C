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
#include "reconstructionSchemes.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(directEvaporation, 0);
    addToRunTimeSelectionTable(macroModel,directEvaporation, components);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::directEvaporation::directEvaporation
(
    const phaseModel& phase1,
    const phaseModel& phase2,
    const volScalarField& p,
    singleComponentSatProp& satModel,
    const compressibleInterPhaseTransportModel& turbModel,
    const dictionary& dict
)
:
    macroModel
    (
        typeName,
        phase1,
        phase2,
        p,
        satModel,
        turbModel,
        dict
    ),
    superheated_(modelDict().get<scalar>("superheated")),
    boilingHeatFlow_
    (
        IOobject
        (
            "boilingHeatFlow",
            p.mesh().time().timeName(),
            p.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        p.mesh(),
        dimensionedScalar("0",dimPower/dimVol,0)
    )
{
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //


// ************************************************************************* //

void Foam::directEvaporation::TSource1(fvScalarMatrix& T1Eqn)
{
    Info << " Tsource " << endl;
    const volScalarField& TSat = satModel_.TSat();
    const volScalarField& T1 = phase1_.thermo().T();

    const fvMesh& mesh = T1.mesh();

    // assume that alpha1 is the more dense phase
    dimensionedScalar superHeated("superHeated",dimTemperature,superheated_);
    dimensionedScalar deltaT = mesh.time().deltaT();

    boilingHeatFlow_ =
        phase1_*pos0(T1.oldTime()-(TSat+superHeated))
        *(T1.oldTime()-(TSat+superHeated))
        *phase1_.thermo().Cp()*phase1_.thermo().rho()/deltaT;


    boilingHeatFlow_.correctBoundaryConditions();

    T1Eqn.source() -= boilingHeatFlow_.internalField()*mesh.V();

}


void Foam::directEvaporation::TSource2(fvScalarMatrix& T2Eqn)
{

}



void Foam::directEvaporation::energySource(volScalarField& Q)
{
    Info << " energySource " << endl;
    const volScalarField& T1 = phase1_.thermo().T();
    const fvMesh& mesh = T1.mesh();

    reconstructionSchemes& surf =
        mesh.lookupObjectRef<reconstructionSchemes>("reconstructionScheme");
    surf.reconstruct(false);

    const volVectorField& normal = surf.normal();
    const volScalarField& TSat = satModel_.TSat();

    dimensionedScalar deltaT = mesh.time().deltaT();
    dimensionedScalar superHeated("superHeated",dimTemperature,superheated_);

    boilingHeatFlow_ =
        phase1_*pos0(T1.oldTime()-(TSat+superHeated))
        *(T1.oldTime()-(TSat+superHeated))
        *phase1_.thermo().Cp()*phase1_.thermo().rho()/deltaT;

    dimensionedScalar evapEnergy = fvc::domainIntegrate(boilingHeatFlow_);

    volScalarField interface
    (
        "interface",
        mag(normal)
    );

    volScalarField scaledEvapSource
    (
        interface*dimensionedScalar("0",dimless,1)
    );

    dimensionedScalar area("area",dimArea,gSum(mag(normal)().internalField()));
    if (area.value() != 0)
    {
        dimensionedScalar scale = evapEnergy/area;
        scaledEvapSource *= scale;
    }

    scaledEvapSource.ref() /= mesh.V();

    Q += scaledEvapSource;

}


void Foam::directEvaporation::energySource1(volScalarField& q1)
{

}


void Foam::directEvaporation::energySource2(volScalarField& q2)
{

}


void Foam::directEvaporation::massSource(volScalarField& rhoSource)
{

}

void Foam::directEvaporation::alphaSource(volScalarField& rhoSource)
{

}
