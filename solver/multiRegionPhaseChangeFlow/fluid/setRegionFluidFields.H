    dynamicFvMesh& mesh = fluidRegions[i];

    volScalarField& p_rgh = p_rghFluid[i];

    volVectorField& U = UFluid[i];

    surfaceScalarField& phi = phiFluid[i];

    // alphaPhi
    surfaceScalarField& alphaPhi10 = alphaPhiFluid[i];
    tmp<surfaceScalarField> talphaPhi1Corr0;

    // thermo lib
    twoPhaseModelThermo& mixture = thermoFluid[i];

    volScalarField& alpha1(mixture.alpha1());
    volScalarField& alpha2(mixture.alpha2());

    Info<< "Reading thermophysical properties\n" << endl;

    volScalarField& p = mixture.p();
    volScalarField& T = mixture.T();
    volScalarField& rho1 = mixture.thermo1().rho();
    const volScalarField& psi = mixture.psi();
    const volScalarField& psi1 = mixture.thermo1().psi();
    volScalarField& rho2 = mixture.thermo2().rho();
    const volScalarField& psi2 = mixture.thermo2().psi();

    volScalarField& TL = mixture.thermo1().T();
    volScalarField& TV = mixture.thermo2().T();

    volScalarField& dpdt = dpdtFluid[i];

    advectionSchemes&  advector = advectionSchemesFluid[i];

    volScalarField& rho = rhoFluid[i];

    volScalarField& rhoCpL = rhoCpLFluid[i];
    volScalarField& rhoCpV = rhoCpVFluid[i];

    dimensionedScalar pMin = pMinFluid[i];

    dimensionedScalar TMin = TMinFluid[i];
    dimensionedScalar TMax = TMaxFluid[i];

    const volScalarField& gh = surfForcesFluid[i].acc();

    surfaceScalarField& rhoPhi = rhoPhiFluid[i];

    surfaceScalarField& rhoCpLPhi = rhoCpLPhiFluid[i];
    surfaceScalarField& rhoCpVPhi = rhoCpVPhiFluid[i];

    volScalarField& dgdt = dgdtFluid[i];

    interfaceProperties& interface = interfaceFluid[i];
    bool frozenFlow = frozenFlowFluid[i];

    compressibleInterPhaseTransportModel& turbulence = turbulenceFluid[i];

    surfaceForces& surfForces = surfForcesFluid[i];
    singleComponentPhaseChange& phaseChange = PhaseChangeModelFluid[i];

    dimensionedScalar& initialMass = initialMassFluid[i];

    volScalarField& K = KFluid[i];

    fv::options& fvOptions = fvOptionsFluid[i];
    IOMRFZoneList& MRF = MRFfluid[i];
