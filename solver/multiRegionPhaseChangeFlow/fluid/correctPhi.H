
    if (mesh.changing())
    {
        volVectorField::Boundary& Ubf = U.boundaryFieldRef();
        surfaceScalarField::Boundary& phibf =
            phi.boundaryFieldRef();

        forAll(Ubf, patchi)
        {
            if (Ubf[patchi].fixesValue())
            {
                Ubf[patchi].initEvaluate();
            }
        }

        forAll(Ubf, patchi)
        {
            if (Ubf[patchi].fixesValue())
            {
                Ubf[patchi].evaluate();

                phibf[patchi] = Ubf[patchi] & mesh.Sf().boundaryField()[patchi];
            }
        }
    }

    // Initialize BCs list for pcorr to zero-gradient
    wordList pcorrTypes
    (
        p.boundaryField().size(),
        zeroGradientFvPatchScalarField::typeName
    );

    // Set BCs of pcorr to fixed-value for patches at which p is fixed
    forAll(p_rgh.boundaryField(), patchi) // shoudl be p_rgh not p
    {
        if (p_rgh.boundaryField()[patchi].fixesValue())
        {
            pcorrTypes[patchi] = fixedValueFvPatchScalarField::typeName;
        }
    }

    volScalarField pcorr
    (
        IOobject
        (
            "pcorr",
            runTime.timeName(),
            mesh
        ),
        mesh,
        dimensionedScalar(p.dimensions(), Zero),
        pcorrTypes
    );

    // temporary comment out with it pressurization simulation are not possible
    // if (pcorr.needReference())
    // {
    //     fvc::makeRelative(phi, U);
    //     adjustPhi(phi, U, pcorr);
    //     fvc::makeAbsolute(phi, U);
    // }

    mesh.setFluxRequired(pcorr.name());
    dimensionedScalar rAUf("rAUf", dimTime/rho.dimensions(), 1);

    for(int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
    {
        // Solve for pcorr such that the divergence of the corrected flux
        // matches the divU provided (from previous iteration, time-step...)
        fvScalarMatrix pcorrEqn
        (
            fvm::laplacian(rAUf, pcorr) == fvc::div(phi) - absdivU
        );

        label refCell = -1;
        if(Pstream::master())
        {
            refCell = 0;
        }
        pcorrEqn.setReference(refCell, 0.0);

        pcorrEqn.solve
        (
            mesh.solver(pcorr.select(nonOrth == nNonOrthCorr))
        );

        if (nonOrth == nNonOrthCorr)
        {
            phi -= pcorrEqn.flux();
        }
    }

    //  volScalarField contErr(fvc::div(phi));

    //  scalar sumLocalContErr = runTime.deltaTValue()*
    //      mag(contErr)().weightedAverage(mesh.V()).value();

    //  scalar globalContErr = runTime.deltaTValue()*
    //      contErr.weightedAverage(mesh.V()).value();
    //  cumulativeContErr += globalContErr;

    //  Info<< "time step continuity errors : sum local = " << sumLocalContErr
    //      << ", global = " << globalContErr
    //      << ", cumulative = " << cumulativeContErr
    //      << endl;
//#include "compressibleContinuityErrors.H"
