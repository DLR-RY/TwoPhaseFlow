#include "advectionSchemes.H"
#include "fvCFD.H"
#include "dynamicFvMesh.H"
#include "isoAdvection.H"
#include "EulerDdtScheme.H"
#include "localEulerDdtScheme.H"
#include "CrankNicolsonDdtScheme.H"
#include "subCycle.H"
#include "immiscibleIncompressibleTwoPhaseMixture.H"

#include "turbulentTransportModel.H"
#include "pimpleControl.H"
#include "fvOptions.H"
#include "CorrectPhi.H"
#include "fvcSmooth.H"
#include "surfaceForces.H"

#include "dropletCloud.H"
#include "phaseCoupling.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Solver for two incompressible, isothermal immiscible fluids"
        " using isoAdvector phase-fraction based interface capturing.\n"
        "With optional mesh motion and mesh topology changes including"
        " adaptive re-meshing.\n"
        "The solver is derived from interFoam"
    );

    Foam::argList::addBoolOption
    (
        "overwrite",
        "Update and overwrite the existing mesh useful for adaptive mesh refinement"
    );

    #include "postProcess.H"

    #include "addCheckCaseOptions.H"
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createDynamicFvMesh.H"
    #include "initContinuityErrs.H"
    #include "createDyMControls.H"
    #include "createFields.H"
    #include "initCorrectPhi.H"
    #include "createUfIfPresent.H"

    const bool overwrite = args.found("overwrite");

    turbulence->validate();

    #include "CourantNo.H"
    #include "setInitialDeltaT.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        #include "readDyMControls.H"
        //#include "porousCourantNo.H"
       // #include "porousAlphaCourantNo.H"
        #include "CourantNo.H"
        #include "alphaCourantNo.H"
        #include "setDeltaT.H"

        ++runTime;

        if(overwrite)
        {
            runTime.setTime(runTime.value() - runTime.deltaTValue(), 1);
            runTime.writeAndEnd();
        }

        Info<< "Time = " << runTime.timeName() << nl << endl;

        // Inject droplets
        coupling.update();

        cloud.storeGlobalPositions();

        // Move droplets
        cloud.move();

        // --- Pressure-velocity PIMPLE corrector loop
        while (pimple.loop())
        {
            if (pimple.firstIter() || moveMeshOuterCorrectors)
            {
                advector->surf().reconstruct();

                mesh.update();
            }

            if(overwrite)
            {
                continue;
            }


            #include "alphaControls.H"
            #include "alphaEqnSubCycle.H"

            mixture.correct();

            surfForces.correct();

            if (pimple.frozenFlow())
            {
                continue;
            }

            #include "UEqn.H"

            // --- Pressure corrector loop
            while (pimple.correct())
            {
                #include "pEqn.H"
            }

            if (pimple.turbCorr())
            {
                turbulence->correct();
            }
        }

        runTime.write();

        runTime.printExecutionTime(Info);
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
