// If there are more than one outer corrector, we use a mixture of old and
// new U and phi for propagating alpha1 in all but the first outer iteration
if (oCorr != 0)
{

       alpha1 = alpha1.oldTime();

}

tmp<surfaceScalarField> talphaPhi1(alphaPhi10);
volScalarField tCpl(mixture.thermo1().Cp());
volScalarField tCpv(mixture.thermo2().Cp());

if (nAlphaSubCycles > 1)
{
    dimensionedScalar totalDeltaT = runTime.deltaT();

    talphaPhi1 = new surfaceScalarField
    (
        IOobject
        (
            "alphaPhi1",
            runTime.timeName(),
            mesh
        ),
        mesh,
        dimensionedScalar(alphaPhi10.dimensions(), Zero)
    );

    surfaceScalarField rhoPhiSum
    (
        IOobject
        (
            "rhoPhiSum",
            runTime.timeName(),
            mesh
        ),
        mesh,
        dimensionedScalar(rhoPhi.dimensions(), Zero)
    );

    surfaceScalarField rhoCpLPhiSum
    (
        IOobject
        (
            "rhoPhiSum",
            runTime.timeName(),
            mesh
        ),
        mesh,
        dimensionedScalar("0", rhoCpLPhi.dimensions(), 0)
    );

    surfaceScalarField rhoCpVPhiSum
    (
        IOobject
        (
            "rhoPhiSum",
            runTime.timeName(),
            mesh
        ),
        mesh,
        dimensionedScalar("0", rhoCpVPhi.dimensions(), 0)
    );

    for
    (
        subCycle<volScalarField> alphaSubCycle(alpha1, nAlphaSubCycles);
        !(++alphaSubCycle).end();
    )
    {
        #include "alphaEqn.H"
        talphaPhi1.ref() += (runTime.deltaT()/totalDeltaT)*alphaPhi10;
        rhoPhiSum += (runTime.deltaT()/totalDeltaT)*rhoPhi;
        rhoCpLPhiSum += (runTime.deltaT()/totalDeltaT)*rhoCpLPhi;
        rhoCpVPhiSum += (runTime.deltaT()/totalDeltaT)*rhoCpVPhi;
    }

    rhoPhi = rhoPhiSum;
    rhoCpLPhi = rhoCpLPhiSum;
    rhoCpVPhi = rhoCpVPhiSum;

}
else
{
    #include "alphaEqn.H"
}


rho = alpha1*rho1 + alpha2*rho2;

rhoCpL = alpha1*rho1*tCpl + alpha2*rho1*tCpl*1e-7;
rhoCpV = alpha2*rho2*tCpv + alpha1*rho2*tCpv*1e-7;

const surfaceScalarField& alphaPhi1 = talphaPhi1();
surfaceScalarField alphaPhi2("alphaPhi2", phi - alphaPhi1);
mixture.phase1().alphaPhi() = alphaPhi10;
mixture.phase2().alphaPhi() = phi - alphaPhi10;

volScalarField::Internal contErr
(
    (
        fvc::ddt(rho) + fvc::div(rhoPhi)
      - (fvOptions(alpha1, mixture.thermo1().rho())&rho1)
      - (fvOptions(alpha2, mixture.thermo2().rho())&rho2)
    )()
);
