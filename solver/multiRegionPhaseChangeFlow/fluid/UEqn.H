    fvVectorMatrix UEqn
    (
        fvm::ddt(rho, U) + fvm::div(rhoPhi, U)
      - fvm::Sp(contErr, U)
      + MRF.DDt(rho, U)
      + turbulence.divDevRhoReff(U)
     ==
        fvOptions(rho, U)
    );

    UEqn.relax();

    fvOptions.constrain(UEqn);

    if (momentumPredictor)
    {
        solve
        (
            UEqn
         ==
            fvc::reconstruct
            (
                (
                    //  mixture.surfaceTensionForce()
                    surfForces.surfaceTensionForce()
                    //- ghf*fvc::snGrad(rho)
                    + surfForces.accelerationForce()
                    - fvc::snGrad(p_rgh)
                ) * mesh.magSf()
            ),
            mesh.solver(U.select(finalIter))
        );

        fvOptions.correct(U);
        K = 0.5*magSqr(U);
    }
