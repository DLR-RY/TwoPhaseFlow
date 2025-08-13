/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016-2017 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2016 OpenCFD Ltd.
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

#include "superHeated.H"
#include "surfaceInterpolate.H"
#include "fvcSnGrad.H"
#include "wallPolyPatch.H"
#include "turbulentFluidThermoModel.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(superHeated, 0);
    addToRunTimeSelectionTable(functionObject, superHeated, dictionary);
}
}


// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

void Foam::functionObjects::superHeated::writeFileHeader(Ostream& os) const
{
    // Add headers to output data
    writeHeader(os, "superHeated ");
    writeCommented(os, "Time");
    writeTabbed(os, "maxSuperHeated");
    os  << endl;
}



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::superHeated::superHeated
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    writeFile(obr_, name, typeName, dict),
    TSat_(),
    T_(),
    alpha_()
{
    volScalarField* superHeatedPtr
    (
        new volScalarField
        (
            IOobject
            (
                type(),
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar("0", dimTemperature, 0)
        )
    );

    mesh_.objectRegistry::store(superHeatedPtr);

    read(dict);

    writeFileHeader(file());
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::superHeated::~superHeated()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::superHeated::read(const dictionary& dict)
{
    fvMeshFunctionObject::read(dict);
    writeFile::read(dict);

    TSat_ = word(dict.lookup("TSat"));
    T_ = word(dict.lookup("T"));
    alpha_ = word(dict.lookup("alpha"));

    return true;
}


bool Foam::functionObjects::superHeated::execute()
{
    volScalarField& superHeated = const_cast<volScalarField&>
    (
        lookupObject<volScalarField>(type())
    );

    volScalarField& TSat = const_cast<volScalarField&>
    (
        lookupObject<volScalarField>(TSat_)
    );

    volScalarField& T = const_cast<volScalarField&>
    (
        lookupObject<volScalarField>(T_)
    );

    volScalarField& alpha = const_cast<volScalarField&>
    (
        lookupObject<volScalarField>(alpha_)
    );

    superHeated = pos(alpha-scalar(1e-4))*T - TSat;


    return true;
}


bool Foam::functionObjects::superHeated::write()
{
    const volScalarField& superHeated = lookupObject<volScalarField>(type());

    superHeated.write();

    scalar maximumValue = gMax(superHeated.internalField());

    Log << type() << " " << name() << " write:" << nl
        << "    writing field " << superHeated.name() << endl;

    superHeated.write();

    if (Pstream::master())
    {
        writeCurrentTime(file());

        file()
            << token::TAB << maximumValue
            << endl;
    }



    return true;
}


// ************************************************************************* //
