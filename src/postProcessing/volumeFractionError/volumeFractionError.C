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

#include "volumeFractionError.H"
#include "addToRunTimeSelectionTable.H"

#include "advectionSchemes.H"
#include "cutCellImpFunc.H"
#include "cutCellIso.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(volumeFractionError, 0);
    addToRunTimeSelectionTable(functionObject, volumeFractionError, dictionary);
}
}


// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

void Foam::functionObjects::volumeFractionError::writeFileHeader(Ostream& os) const
{
    // Add headers to output data
    writeHeader(os, "volumeFractionError ");
    writeCommented(os, "Time");
    writeTabbed(os, "Eshape");
    writeTabbed(os, "Emass");
    writeTabbed(os, "Ebound");
    writeTabbed(os, "tRec");
    writeTabbed(os, "Tadv");
    writeTabbed(os, "celltoCellDist");
    writeTabbed(os, "Ncells");
    writeTabbed(os, "avg(CellVol)");
    os  << endl;
}



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::volumeFractionError::volumeFractionError
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    writeFile(obr_, name, typeName, dict),
    Field_(),
    dict_(dict),
    initMass_(0),
    initCentre_(dict.get<vector>("origin"))
{



    read(dict);

    writeFileHeader(file());
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::volumeFractionError::~volumeFractionError()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::volumeFractionError::read(const dictionary& dict)
{
    fvMeshFunctionObject::read(dict);
    writeFile::read(dict);

    Field_ = dict.get<word>("field");


    // set initial mass
    Foam::autoPtr<Foam::implicitFunction> func
    (
        implicitFunction::New
        (
            dict.get<word>("functionType"),
            dict
        )
    );

    scalarField f(mesh_.nPoints(),0.0);


    forAll(f,pI)
    {
        f[pI] =  func->value(mesh_.points()[pI]);
    }


    cutCellIso cutCell(mesh_,f);

    volScalarField alphaExact
    (
          IOobject
          (
              "alphaExact",
              mesh_.time().timeName(),
              mesh_,
              IOobject::NO_READ,
              IOobject::NO_WRITE
          ),
          mesh_,
          dimensionedScalar("0", dimless, 0),
          "calculated"
    );


    forAll(alphaExact,cellI)
    {
        label cellStatus = cutCell.calcSubCell(cellI,0.0);

        if(cellStatus == -1)
        {
            alphaExact[cellI] = 1;
        }
        else if(cellStatus == 1)
        {
            alphaExact[cellI] = 0;
        }
        else if(cellStatus == 0)
        {
            alphaExact[cellI]= cutCell.VolumeOfFluid();
        }

    }
    const volScalarField& alpha = lookupObject<volScalarField>(Field_);

    initMass_ = fvc::domainIntegrate(alpha.internalField()).value();

    return true;
}


bool Foam::functionObjects::volumeFractionError::execute()
{
    // do nothing

    return true;
}


bool Foam::functionObjects::volumeFractionError::write()
{

    const volScalarField& alpha = lookupObject<volScalarField>(Field_);

    vector centre = dict_.get<vector>("origin");

    const volVectorField& U = mesh_.lookupObject<volVectorField>("U");

    volScalarField magU ("magU",mag(U));

    scalar maxMagU = max(magU).value();
    scalar minMagU = min(magU).value();

    reduce(maxMagU,maxOp<scalar>());
    reduce(minMagU,minOp<scalar>());

    vector U0 = Foam::vector(0,0,0);

    if((maxMagU - minMagU) <= SMALL )
    {
        U0 = U[0];

    }

    vector centrePos = initCentre_ + U0*mesh_.time().value();






    dict_.set<vector>("origin",centrePos);

    Foam::autoPtr<Foam::implicitFunction> func =  implicitFunction::New
    (
           word(dict_.lookup("functionType")),
           dict_
    );



    scalarField f(mesh_.nPoints(),0.0);

    forAll(f,pI)
    {
        f[pI] =  func->value(mesh_.points()[pI]);
    }


     cutCellIso cutCell(mesh_,f); // changes f

    volScalarField alphaExact
    (
          IOobject
          (
              "alphaExact",
              mesh_.time().timeName(),
              mesh_,
              IOobject::NO_READ,
              IOobject::NO_WRITE
          ),
          mesh_,
          dimensionedScalar("0", dimless, 0),
          "calculated"
    );

//    cutCellImpFunc cutCell(mesh_,f,func); // changes f


    forAll(alphaExact,cellI)
    {
        label cellStatus = cutCell.calcSubCell(cellI,0.0);

        if(cellStatus == -1)
        {
            alphaExact[cellI] = 1;
        }
        else if(cellStatus == 1)
        {
            alphaExact[cellI] = 0;
        }
        else if(cellStatus == 0)
        {
            alphaExact[cellI]= cutCell.VolumeOfFluid();
        }

    }


    scalar E1 = fvc::domainIntegrate(mag(alpha.internalField() - alphaExact.internalField())).value();
    scalar Emass = mag(fvc::domainIntegrate(alpha.internalField()).value() - initMass_);

    scalarField fluidVol = alpha.primitiveField();
    scalar Ebound = mag(max(-min(fluidVol),max(fluidVol-1)));

    reduce(Ebound,maxOp<scalar>());





    const advectionSchemes& advection = mesh_.lookupObjectRef<advectionSchemes>("advectionSchemes");

    scalar reconAvgTimePerIter = advection.runTime().x() / mesh_.time().timeIndex();
    scalar advectionAvgTimePerIter = advection.runTime().y() / mesh_.time().timeIndex();

    scalar length = average(mag(mesh_.delta())()).value();
    scalar nCells = mesh_.nCells();
    scalar avgVol = average(mesh_.V()).value();


//    Log << type() << " " << volumeFractionError.name() << " is " << intValue << endl;


    if (Pstream::master())
    {
        writeCurrentTime(file());

        file()
            << token::TAB << E1
            << token::TAB << Emass
            << token::TAB << Ebound
            << token::TAB << reconAvgTimePerIter
            << token::TAB << advectionAvgTimePerIter
            << token::TAB << length
            << token::TAB << nCells
            << token::TAB << avgVol
            << endl;
    }



    return true;
}


// ************************************************************************* //
