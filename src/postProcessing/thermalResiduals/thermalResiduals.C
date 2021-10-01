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

#include "thermalResiduals.H"
#include "surfaceInterpolate.H"
#include "fvcSnGrad.H"
#include "wallPolyPatch.H"
#include "turbulentFluidThermoModel.H"
#include "addToRunTimeSelectionTable.H"
#include "singleComponentPhaseChange.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(thermalResiduals, 0);
    addToRunTimeSelectionTable(functionObject, thermalResiduals, dictionary);
}
}


// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

void Foam::functionObjects::thermalResiduals::writeFileHeader(Ostream& os) const
{
    // Add headers to output data
    writeHeader(os, "thermalResiduals ");
    writeCommented(os, "Time");
    writeTabbed(os, "thermalResiduals_Vol");
    writeTabbed(os, "thermalResiduals_Mass");
    os  << endl;
}



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::thermalResiduals::thermalResiduals
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    writeFile(obr_, name, typeName, dict),
    TField_(),
    alphaField_(),
    rhoField_(),
    subCooled_("subCooled",dimTemperature,dict.get<scalar>("subCooled"))
{



    read(dict);

    writeFileHeader(file());
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::thermalResiduals::~thermalResiduals()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::thermalResiduals::read(const dictionary& dict)
{
    fvMeshFunctionObject::read(dict);
    writeFile::read(dict);

    TField_ = word(dict.lookup("T"));
    alphaField_ = word(dict.lookup("alpha"));
    rhoField_ = word(dict.lookup("rho"));


    return true;
}


bool Foam::functionObjects::thermalResiduals::execute()
{
    // do nothing

    return true;
}


bool Foam::functionObjects::thermalResiduals::write()
{

    const volScalarField& T = lookupObject<volScalarField>(TField_);
    singleComponentPhaseChange& pcModel = lookupObjectRef<singleComponentPhaseChange>("singleComponentPhaseChange");
    const volScalarField& alpha1 = lookupObject<volScalarField>(alphaField_);
    const volScalarField& rho1 = lookupObject<volScalarField>(rhoField_);
    const volScalarField& p = lookupObject<volScalarField>("p");
    const volScalarField& TSat = pcModel.satProp().TSat();


    volScalarField thermalResVol = neg0(TSat-(T+subCooled_))*alpha1;
    volScalarField thermalResMass = neg0(TSat-(T+subCooled_))*alpha1*rho1;


    scalar intVolValue = fvc::domainIntegrate(thermalResVol.internalField()).value();
    scalar intMassValue = fvc::domainIntegrate(thermalResMass.internalField()).value();


    //Log << type() << " " << thermalResiduals.name() << " are " << intValue << endl;


    if (Pstream::master())
    {
        writeCurrentTime(file());

        file()
            << token::TAB << intVolValue
            << token::TAB << intMassValue
            << endl;
    }



    return true;
}


// ************************************************************************* //
