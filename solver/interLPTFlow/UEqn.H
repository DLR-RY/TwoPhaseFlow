    MRF.correctBoundaryVelocity(U);

    fvVectorMatrix UEqn
    (
        fvm::ddt(rho, U) + fvm::div(rhoPhi, U)
      + MRF.DDt(rho, U)
      + turbulence->divDevRhoReff(rho, U)
      + fvm::Sp(rho*coupling.damping(), U)
      + fvm::Su(rho*coupling.source(), U)
     ==
        fvOptions(rho, U)
    );

    #include "UEqnAddPorosity.H"

    UEqn.relax();

    fvOptions.constrain(UEqn);

    if (pimple.momentumPredictor())
    {
        solve
        (
            UEqn
         ==
            fvc::reconstruct
            (
                (
                    //mixture.surfaceTensionForce()
                  surfForces.surfaceTensionForce()
                  //- ghf*fvc::snGrad(rho)
                  + surfForces.accelerationForce()
                  - fvc::snGrad(p_rgh)
                ) * mesh.magSf()
            )
        );

        fvOptions.correct(U);
    }
