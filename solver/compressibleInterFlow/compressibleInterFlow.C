/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2017 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2017 OpenCFD Ltd
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

Application
    compressibleInterFlow

Description
    Solver for two compressible, non-isothermal immiscible fluids using a VOF
    (volume of fluid) phase-fraction based interface capturing approach,
    with optional mesh motion and mesh topology changes including adaptive
    re-meshing.

    The momentum and other fluid properties are of the "mixture" and a single
    momentum equation is solved.

    Turbulence modelling is generic, i.e. laminar, RAS or LES may be selected.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "dynamicFvMesh.H"
#include "CMULES.H"
#include "EulerDdtScheme.H"
#include "localEulerDdtScheme.H"
#include "CrankNicolsonDdtScheme.H"
#include "subCycle.H"
#include "compressibleInterPhaseTransportModel.H"
#include "pimpleControl.H"
#include "fvOptions.H"
#include "CorrectPhi.H"
#include "fvcSmooth.H"
#include "compressibleInterPhaseTransportModel.H"

#include "advectionSchemes.H"
#include "surfaceForces.H"
#include "twoPhaseModelThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Solver for two compressible, non-isothermal immiscible fluids"
        " using VOF phase-fraction based interface capturing.\n"
        "With optional mesh motion and mesh topology changes including"
        " adaptive re-meshing."
    );

    Foam::argList::addBoolOption
    (
        "overwrite",
        "Update and overwrite the existing mesh useful for adaptive mesh refinement"
    );

    #include "postProcess.H"

    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createDynamicFvMesh.H"
    #include "initContinuityErrs.H"
    #include "createDyMControls.H"
    #include "createFields.H"
    #include "createUf.H"
    #include "CourantNo.H"
    #include "setInitialDeltaT.H"

    const bool overwrite = args.found("overwrite");

    volScalarField& p = mixture.p();
    volScalarField& T = mixture.T();
    const volScalarField& psi1 = mixture.thermo1().psi();
    const volScalarField& psi2 = mixture.thermo2().psi();

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        #include "readDyMControls.H"

        // Store divU from the previous mesh so that it can be mapped
        // and used in correctPhi to ensure the corrected phi has the
        // same divergence
        volScalarField divU("divU0", fvc::div(fvc::absolute(phi, U)));

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

        // --- Pressure-velocity PIMPLE corrector loop
        while (pimple.loop())
        {
            if (pimple.firstIter() || moveMeshOuterCorrectors)
            {
                scalar timeBeforeMeshUpdate = runTime.elapsedCpuTime();

                advector->surf().reconstruct();

                mesh.update();

                if (mesh.changing())
                {
                    MRF.update();

                    // advector->surf().mapAlphaField();
                    // alpha2 = 1.0 - alpha1;
                    // alpha2.correctBoundaryConditions();
                    // rho == alpha1*rho1 + alpha2*rho2;
                    // rho.correctBoundaryConditions();
                    // rho.oldTime() = rho;
                    // alpha2.oldTime() = alpha2;

                    Info<< "Execution time for mesh.update() = "
                        << runTime.elapsedCpuTime() - timeBeforeMeshUpdate
                        << " s" << endl;

                    // gets recompute by surfaces forces
                    // gh = (g & mesh.C()) - ghRef;
                    // ghf = (g & mesh.Cf()) - ghRef;
                }

                if ((mesh.changing() && correctPhi))
                {
                    // Calculate absolute flux from the mapped surface velocity
                    phi = mesh.Sf() & Uf;

                    #include "correctPhi.H"

                    // Make the fluxes relative to the mesh motion
                    fvc::makeRelative(phi, U);

                    mixture.correct();
                }

                if (mesh.changing() && checkMeshCourantNo)
                {
                    #include "meshCourantNo.H"
                }
            }

            if(overwrite)
            {
                continue;
            }

            #include "alphaControls.H"
            #include "compressibleAlphaEqnSubCycle.H"

            mixture.correct();

            turbulence.correctPhasePhi();
            surfForces.correct();
            #include "YEqns.H"

            #include "UEqn.H"
            #include "TEqn.H"

            // --- Pressure corrector loop
            while (pimple.correct())
            {
                #include "pEqn.H"
            }

            if (pimple.turbCorr())
            {
                turbulence.correct();
            }
        }

        rho = alpha1*rho1 + alpha2*rho2;

        // Correct p_rgh for consistency with p and the updated densities
        p_rgh = p - rho*gh;
        p_rgh.correctBoundaryConditions();

        runTime.write();

        runTime.printExecutionTime(Info);
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
