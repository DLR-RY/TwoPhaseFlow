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

#include "interfaceRegion.H"
#include "reconstructionSchemes.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(interfaceRegion, 0);
    addToRunTimeSelectionTable(functionObject, interfaceRegion, dictionary);
}
}




// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::interfaceRegion::interfaceRegion
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    markIR_(mesh_),
    interfaceRegion_
    (
        IOobject
        (
            "interfaceRegion",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("0", dimless, 0)
    ),
    nLayers_(0)
{

    read(dict);

    execute();

}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::interfaceRegion::~interfaceRegion()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::interfaceRegion::read(const dictionary& dict)
{
    fvMeshFunctionObject::read(dict);

    nLayers_ = dict.get<scalar>("nLayers");


    return true;
}


bool Foam::functionObjects::interfaceRegion::execute()
{
    // do nothing
    reconstructionSchemes& surf = mesh_.lookupObjectRef<reconstructionSchemes>("reconstructionScheme");

    boolList nextToInterface(mesh_.nCells(),false);
    markIR_.markCellsNearSurf
    (
        surf.interfaceCell(),
        nLayers_,
        nextToInterface,
        interfaceRegion_
    );



    return true;
}


bool Foam::functionObjects::interfaceRegion::write()
{
    return true;
}


// ************************************************************************* //
