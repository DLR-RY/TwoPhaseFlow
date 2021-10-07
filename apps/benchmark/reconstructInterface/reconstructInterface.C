/*---------------------------------------------------------------------------*\
            Copyright (c) 2017-2021, German Aerospace Center (DLR)
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
    initAlphaField

Description
    Uses cellCellIso to create a volume fraction field from either a cylinder,
    a sphere or a plane.

Author
    Henning Scheufler, DLR, all rights reserved.
    Johan Roenby, DHI, all rights reserved.


\*---------------------------------------------------------------------------*/

#include "fvCFD.H"

#include "OFstream.H"

#include "reconstructionSchemes.H"
#include "implicitFunction.H"
#include "cutCellImpFunc.H"
#include "cutCellIso.H"
#include "reconstructionError.H"

#include "surfaceForces.H"




// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class T>
void setAlpha
(
    T& cutCell,
    const fvMesh& mesh,
    const dictionary& initAlphaFieldDict,
    volScalarField& alpha1
)
{

    forAll(alpha1,cellI)
    {
        label cellStatus = cutCell.calcSubCell(cellI,0.0);

        if(cellStatus == -1)
        {
            alpha1[cellI] = 1;
        }
        else if(cellStatus == 1)
        {
            alpha1[cellI] = 0;
        }
        else if(cellStatus == 0)
        {
            if(cellStatus == 0 && mag(cutCell.faceArea()) != 0)
            {
                alpha1[cellI]= max(min(cutCell.VolumeOfFluid(),1),0);
            }
        }
    }

    alpha1.correctBoundaryConditions();

}

int main(int argc, char *argv[])
{
//    #include "addRegionOption.H"
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createControl.H"


    IOdictionary initAlphaFieldDict
    (
        IOobject
        (
            "initAlphaFieldDict",
            mesh.time().system(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );



// init
    Info<< "Reading field U\n" << endl;
    volVectorField U
    (
        IOobject
        (
            "U",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        mesh
    );

    #include "createPhi.H"

    Info<< "Reading field alpha1\n" << endl;
    volScalarField alpha1
    (
        IOobject
        (
            "alpha1",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    IOdictionary transportProperties
    (
        IOobject
        (
            "transportProperties",
            U.time().constant(),
            U.db(),
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    );

    reconstructionError recErr( mesh,mesh,initAlphaFieldDict);

    #include "createTimeControls.H"
    #include "readTimeControls.H"
    #include "CourantNo.H"
    #include "setInitialDeltaT.H"

    //Setting velocity field and face fluxes for next time step
    IOdictionary fvSolutionDict
    (
        IOobject
        (
            "fvSolution",
            alpha1.time().system(),
            alpha1.db(),
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    );

    vector centre = initAlphaFieldDict.get<vector>("origin");
    //    scalar radius = readScalar(initAlphaFieldDict.lookup("radius"));

    Random rndCentre(1234567);

    autoPtr<reconstructionSchemes> surf =
        reconstructionSchemes::New(alpha1,phi,U,fvSolutionDict);
    surfaceForces surfForces(alpha1,phi,U,transportProperties);

    //    label nIter = 100;


    dictionary reconDict = fvSolutionDict.subDict("reconstruction");

    // lookup of relevevant parameters
    label nIter = reconDict.get<label>("nIter");
    word setAlphaMethod = reconDict.get<word>("setAlphaMethod");
    if (setAlphaMethod != "cutCellImpFunc" && setAlphaMethod != "cutCellIso")
    {
        FatalError  << "valid choice are only cutCellImpFunc or cutCellIso"
                    << abort(FatalError);
    }

    word functionType (initAlphaFieldDict.lookup("type"));
    bool twoDim = (functionType == "cylinder");
    Info << "twoDim = " << twoDim << endl;
    Info << "functionType = " << functionType << endl;

    scalar recTime = 0;
    vector centreMin = centre - 0.1*centre;
    vector centreMax = centre + 0.1*centre;



    while (runTime.run())
    {

        runTime++;

        for(int iteration = 0;iteration < nIter;iteration++)
        {

            vector centrePos = rndCentre.globalPosition < vector > (centreMin,centreMax); //vector::zero;


            if (nIter > 1)
            {
                initAlphaFieldDict.set<vector>("origin",centrePos);
            }

            centre = initAlphaFieldDict.get<vector>("origin");

            Info << "centre " <<  centre << endl;

            Foam::autoPtr<Foam::implicitFunction> func
            (
                implicitFunction::New
                (
                    word(initAlphaFieldDict.get<word>("type")),
                    initAlphaFieldDict
                )
            );

            scalarField f(mesh.nPoints(),0.0);

            forAll(f,pI)
            {
                f[pI] =  func->value(mesh.points()[pI]);
            }


            if (setAlphaMethod == "cutCellImpFunc")
            {
                cutCellImpFunc cutCell(mesh,f,func.ref());
                setAlpha(cutCell,mesh,initAlphaFieldDict,alpha1);
            }
            else
            {
                cutCellIso cutCell(mesh,f);
                setAlpha(cutCell,mesh,initAlphaFieldDict,alpha1);
            }

            mesh.time().cpuTimeIncrement();


            surf->reconstruct();

            surfForces.correct();

            surfForces.surfaceTensionForce(); // update curvature

            recTime += mesh.time().cpuTimeIncrement();
            volScalarField curv = mesh.lookupObjectRef<volScalarField>("K_");


            recErr.calcError(initAlphaFieldDict,surf->centre(),surf->normal(),curv,false);

        }

        recErr.write(recTime/nIter,scalar(0));


        runTime.write();

    }


    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
