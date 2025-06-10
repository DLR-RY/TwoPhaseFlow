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

#include "interfaceEnergyFluxes.H"
#include "surfaceInterpolate.H"
#include "fvcSnGrad.H"
#include "wallPolyPatch.H"
#include "turbulentFluidThermoModel.H"
#include "addToRunTimeSelectionTable.H"
#include "reconstructionSchemes.H"
#include "singleComponentPhaseChange.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(interfaceEnergyFluxes, 0);
    addToRunTimeSelectionTable(functionObject, interfaceEnergyFluxes, dictionary);
}
}


// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

void Foam::functionObjects::interfaceEnergyFluxes::writeFileHeader(Ostream& os) const
{
    // Add headers to output data
    writeHeader(os, "interfaceEnergyFluxes ");
    writeCommented(os, "Time");
    writeTabbed(os, "energyFlux_Liquid");
    writeTabbed(os, "energyFlux_Gas");
    os  << endl;
}


Foam::scalar Foam::functionObjects::interfaceEnergyFluxes::computeTotalEnergy
(
    const volScalarField& energyFlux,
    const volVectorField& normal
)
{
    scalar intEnergy = 0;
    forAll(energyFlux,celli)
    {
        intEnergy += energyFlux[celli]*mag(normal[celli]);
    }

    reduce(intEnergy,sumOp<scalar>());

    return intEnergy;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::interfaceEnergyFluxes::interfaceEnergyFluxes
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    writeFile(obr_, name, typeName, dict),
    interfaceEnergyFluxLiquid_
    (
        IOobject
        (
            "interfaceEnergyFluxLiquid",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("0", dimPower/dimArea, 0)
    ),
    interfaceEnergyFluxGas_
    (
        IOobject
        (
            "interfaceEnergyFluxGas",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("0", dimPower/dimArea, 0)
    )
{
    read(dict);

    writeFileHeader(file());
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::interfaceEnergyFluxes::~interfaceEnergyFluxes()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::interfaceEnergyFluxes::read(const dictionary& dict)
{
    fvMeshFunctionObject::read(dict);
    writeFile::read(dict);

    return true;
}


bool Foam::functionObjects::interfaceEnergyFluxes::execute()
{

    return true;
}


bool Foam::functionObjects::interfaceEnergyFluxes::write()
{
    singleComponentPhaseChange& PCModel =
    lookupObjectRef<singleComponentPhaseChange>("singleComponentPhaseChange");

    reconstructionSchemes& surf =
    lookupObjectRef<reconstructionSchemes>("reconstructionScheme");
    const volVectorField& normal = surf.normal();

    interfaceEnergyFluxLiquid_ = PCModel.energyFlux1();
    interfaceEnergyFluxGas_ = PCModel.energyFlux2();

    if(interfaceEnergyFluxLiquid_.mesh().time().writeTime())
    {
        interfaceEnergyFluxLiquid_.write();
        interfaceEnergyFluxGas_.write();
    }

    scalar intEnergyLiquid = computeTotalEnergy(interfaceEnergyFluxLiquid_,normal);
    scalar intEnergyGas = computeTotalEnergy(interfaceEnergyFluxGas_,normal);

    if (Pstream::master())
    {
        writeCurrentTime(file());

        file()
            << token::TAB << intEnergyLiquid
            << token::TAB << intEnergyGas
            << endl;
    }



    return true;
}


// ************************************************************************* //
