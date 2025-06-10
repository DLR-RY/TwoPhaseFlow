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

Application
    setFieldfromTable

Description
    set field from tabulated data.

\*---------------------------------------------------------------------------*/

#include "List.H"
#include "vector.H"
#include "Vector2D.H"
#include "Tuple2.H"
#include "OFstream.H"
#include "IFstream.H"
#include "csvTableReader.H"
#include "interpolationTable.H"

#include "volFields.H"

#include "fvMesh.H"

#include "argList.H"

#include "argList.H"
#include "Time.H"
#include "implicitFunction.H"



using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    #include "addRegionOption.H"
    #include "setRootCase.H"
    #include "createTime.H"
//  #include "createMesh.H"
    #include "createNamedMesh.H"

    argList::noParallel();

    IOdictionary setFieldTableDict
    (
        IOobject
        (
            "setFieldfromTableDict",
            runTime.system(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );


    const word nameField (setFieldTableDict.lookup("field"));
    const bool writeDistanceField = setFieldTableDict.lookupOrDefault("writeDistanceField",false);
    const bool absolute = setFieldTableDict.lookupOrDefault("absolute",false);
    const scalar minDistValue = setFieldTableDict.lookupOrDefault("minDist",-GREAT);
    const scalar maxDistValue = setFieldTableDict.lookupOrDefault("maxDist",GREAT);

    Info<< "Reading field " <<  nameField << endl;
    Info<< "maxDistValue field " <<  maxDistValue << endl;
    Info<< "minDistValue field " <<  minDistValue << endl;

    volScalarField Field
    (
        IOobject
        (
            nameField,
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ
        ),
        mesh
    );


    Info<< "Reading distance function" << endl;

    Foam::autoPtr<Foam::implicitFunction> func =  implicitFunction::New
    (
           setFieldTableDict.get<word>("type"),
           setFieldTableDict
    );


    volScalarField distanceField
    (
        IOobject
        (
            "distanceField",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("0",dimLength,0.0),
        "calculated"
    );


    forAll(distanceField,celli)
    {
        distanceField[celli] = func->value(mesh.C()[celli]);
    };


    forAll(distanceField.boundaryField(), patchI)
    {
        forAll(distanceField.boundaryField()[patchI],i)
        {
            if(distanceField.boundaryField()[patchI].fixesValue())
            {
                vector c = mesh.C().boundaryField()[patchI][i];

                distanceField.boundaryFieldRef()[patchI][i] = func->value(c);
            }
        }
    }

    if(absolute)
    {
        distanceField = mag(distanceField);
    }

    interpolationTable<scalar> interpolationData(setFieldTableDict);
    forAll(Field, cellI)
    {

        scalar dist = distanceField[cellI];
        if(dist > minDistValue  && dist  < maxDistValue)
        {
            Field[cellI]=interpolationData(dist);
        }

    }

    forAll(Field.boundaryField(), patchI)
    {
        forAll(Field.boundaryField()[patchI],i)
        {
            if(Field.boundaryField()[patchI].fixesValue())
            {
                scalar dist = distanceField.boundaryFieldRef()[patchI][i];
                if(minDistValue < dist  && dist  < maxDistValue)
                {
                    Field.boundaryFieldRef()[patchI][i] = interpolationData(dist);
                }
            }
        }
    }

    Info<< "write field " <<  nameField << endl;
    Field.write();

    if(writeDistanceField)
    {
        distanceField.write();
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
