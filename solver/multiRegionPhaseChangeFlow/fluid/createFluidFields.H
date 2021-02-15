// Initialise fluid field pointer lists
// identisch mit compressiblePhaseChangeFoamHardtWondra createField

PtrList<volScalarField> p_rghFluid(fluidRegions.size());
PtrList<volVectorField> UFluid(fluidRegions.size());
PtrList<surfaceScalarField> phiFluid(fluidRegions.size());
List<autoPtr<surfaceVectorField>> UfFluid(fluidRegions.size());

PtrList<twoPhaseModelThermo> thermoFluid(fluidRegions.size()); // add

PtrList<volScalarField> pFluid(fluidRegions.size());
PtrList<volScalarField> rhoFluid(fluidRegions.size()); // add

PtrList<volScalarField> dpdtFluid(fluidRegions.size()); // add

// PhaseChange and Isoadvector
PtrList<advectionSchemes> advectionSchemesFluid(fluidRegions.size());

// rhoCp
PtrList<volScalarField> rhoCpLFluid(fluidRegions.size());
PtrList<volScalarField> rhoCpVFluid(fluidRegions.size());

// surface ScalarField

PtrList<surfaceScalarField> rhoPhiFluid(fluidRegions.size());
PtrList<surfaceScalarField> rhoCpLPhiFluid(fluidRegions.size());
PtrList<surfaceScalarField> rhoCpVPhiFluid(fluidRegions.size());
PtrList<surfaceScalarField> alphaPhiFluid(fluidRegions.size());

// gravitation and turb and misc

PtrList<uniformDimensionedVectorField> gFluid(fluidRegions.size());
PtrList<compressibleInterPhaseTransportModel> turbulenceFluid(fluidRegions.size());
PtrList<volScalarField> ghFluid(fluidRegions.size());
PtrList<surfaceScalarField> ghfFluid(fluidRegions.size());
PtrList<volScalarField> dgdtFluid(fluidRegions.size());
PtrList<volScalarField> KFluid(fluidRegions.size());
PtrList<interfaceProperties> interfaceFluid(fluidRegions.size());

// Phase Change
PtrList<singleComponentPhaseChange> PhaseChangeModelFluid(fluidRegions.size());
PtrList<IOdictionary> thermophysicalDictFluid(fluidRegions.size());
PtrList<surfaceForces> surfForcesFluid(fluidRegions.size());

//scalars

PtrList<dimensionedScalar> pMinFluid(fluidRegions.size());
PtrList<dimensionedScalar> TMinFluid(fluidRegions.size());
PtrList<dimensionedScalar> TMaxFluid(fluidRegions.size());

PtrList<dimensionedScalar> initialMassFluid(fluidRegions.size());
List<bool> frozenFlowFluid(fluidRegions.size(), false);

PtrList<fv::options> fvOptionsFluid(fluidRegions.size());
PtrList<IOMRFZoneList> MRFfluid(fluidRegions.size());

// Populate fluid field pointer lists
forAll(fluidRegions, i)
{
    Info<< "*** Reading fluid mesh thermophysical properties for region "
        << fluidRegions[i].name() << nl << endl;
    p_rghFluid.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "p_rgh",
                runTime.timeName(),
                fluidRegions[i],
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            fluidRegions[i]
        )
    );

    Info<< "    Adding to UFluid\n" << endl;
    UFluid.set
    (
        i,
        new volVectorField
        (
            IOobject
            (
                "U",
                runTime.timeName(),
                fluidRegions[i],
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            fluidRegions[i]
        )
    );

    Info<< "    Adding to phiFluid\n" << endl;
    phiFluid.set
    (
        i,
        new surfaceScalarField
        (
            IOobject
            (
                "phi",
                runTime.timeName(),
                fluidRegions[i],
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            linearInterpolate(UFluid[i])
                & fluidRegions[i].Sf()
        )
    );

    forAll(UfFluid,i)
    {
        if (fluidRegions[i].dynamic())
        {
            Info<< "Constructing face velocity Uf\n" << endl;

            UfFluid[i].reset
            (
                new surfaceVectorField
                (
                    IOobject
                    (
                        "Uf",
                        runTime.timeName(),
                        fluidRegions[i],
                        IOobject::READ_IF_PRESENT,
                        IOobject::AUTO_WRITE
                    ),
                    fvc::interpolate(UFluid[i])
                )
            );
        }
    }

    thermoFluid.set
    (
        i,
        new twoPhaseModelThermo(UFluid[i],phiFluid[i])
    );

    dpdtFluid.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "dpdt",
                runTime.timeName(),
                fluidRegions[i],
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            fluidRegions[i],
            dimensionedScalar("dpdt",thermoFluid[i].p().dimensions()/dimTime, 0)
        )
    );

    Info<< "    Adding to rhoFluid\n" << endl;
    rhoFluid.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "rho",
                runTime.timeName(),
                fluidRegions[i],
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            thermoFluid[i].alpha1()*thermoFluid[i].thermo1().rho()
            + thermoFluid[i].alpha2()*thermoFluid[i].thermo2().rho()
        )
    );
    rhoCpLFluid.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "rhoCpL",
                runTime.timeName(),
                fluidRegions[i],
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            thermoFluid[i].alpha1()*thermoFluid[i].thermo1().rho()
            *thermoFluid[i].thermo1().Cp()  + thermoFluid[i].alpha2()
            *thermoFluid[i].thermo1().rho()*thermoFluid[i].thermo1().Cp()*1e-7
        )
    );
    rhoCpVFluid.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "rhoCpV",
                runTime.timeName(),
                fluidRegions[i],
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            thermoFluid[i].alpha2()*thermoFluid[i].thermo2().rho()
            *thermoFluid[i].thermo2().Cp() + thermoFluid[i].alpha1()
            *thermoFluid[i].thermo2().rho()*thermoFluid[i].thermo2().Cp()*1e-7
        )
    );

    pMinFluid.set
    (
        i,
        new dimensionedScalar(thermoFluid[i].lookup("pMin"))
    );

    TMaxFluid.set
    (
        i,
        new dimensionedScalar(thermoFluid[i].lookup("TMax"))
    );

    TMinFluid.set
    (
        i,
        new dimensionedScalar(thermoFluid[i].lookup("TMin"))
    );

    Info<< "    creating surfaceForcesModel \n" << endl;
    thermophysicalDictFluid.set
    (
        i,
        new IOdictionary
        (
            IOobject
            (
                "thermophysicalProperties",
                UFluid[i].time().constant(),
                UFluid[i].db(),
                IOobject::MUST_READ_IF_MODIFIED,
                IOobject::NO_WRITE
            )
        )

    );

    surfForcesFluid.set
    (
        i,
        new surfaceForces
        (
            thermoFluid[i].alpha1(),
            phiFluid[i],
            UFluid[i],
            thermophysicalDictFluid[i]
        )
    );

    rhoPhiFluid.set
    (
        i,
        new surfaceScalarField
        (
            IOobject
            (
                "rhoPhi",
                runTime.timeName(),
                fluidRegions[i],
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            fvc::interpolate(rhoFluid[i])*phiFluid[i]
        )
    );

    rhoCpLPhiFluid.set
    (
        i,
        new surfaceScalarField
        (
            IOobject
            (
                "rhoCpLPhi",
                runTime.timeName(),
                fluidRegions[i],
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            fvc::interpolate(rhoFluid[i]*thermoFluid[i].thermo1().Cp())
            *phiFluid[i] // init
        )
    );

    rhoCpVPhiFluid.set
    (
        i,
        new surfaceScalarField
        (
            IOobject
            (
                "rhoCpVPhi",
                runTime.timeName(),
                fluidRegions[i],
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            fvc::interpolate(rhoFluid[i]*thermoFluid[i].thermo1().Cp())
            *phiFluid[i] // init
        )
    );

    // set dgdt to zero may be a better option
    // would be to save it for the restart
    dgdtFluid.set
    (
        i,
        new volScalarField
        (
            "dgdt",
            thermoFluid[i].alpha1()*fvc::div(phiFluid[i])*0
        )
    );

    interfaceFluid.set
    (
        i,
        new interfaceProperties
        (
            thermoFluid[i].alpha1(),
            UFluid[i],
            thermoFluid[i]
        )
    );

    // createAlphaFluxes
    IOobject alphaPhi10Header
    (
        "alphaPhi10",
        runTime.timeName(),
        fluidRegions[i],
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    );

    const bool alphaRestart =
    alphaPhi10Header.typeHeaderOk<surfaceScalarField>(true);

    // MULES flux from previous time-step
    alphaPhiFluid.set
    (
        i,
        new surfaceScalarField
        (
            alphaPhi10Header,
            phiFluid[i]*fvc::interpolate(thermoFluid[i].alpha1())
        )
    );
    // // limit clip alpha1
    // thermoFluid[i].alpha1() =
    //     thermoFluid[i].alpha1()
    //     *pos0(thermoFluid[i].alpha1() - 1e-4)
    //     *neg0(thermoFluid[i].alpha1() - (1.0 - 1e-4))
    //     + pos0(thermoFluid[i].alpha1() - (1.0 - 1e-4));

    advectionSchemesFluid.set
    (
        i,
        advectionSchemes::New(thermoFluid[i].alpha1(),phiFluid[i],UFluid[i])
    );

    Info<< "    Adding turbulence \n" << endl;
    turbulenceFluid.set
    (
        i,
        new compressibleInterPhaseTransportModel
        (
            rhoFluid[i],
            UFluid[i],
            phiFluid[i],
            rhoPhiFluid[i],
            alphaPhiFluid[i],
            thermoFluid[i]
        )
    );

    initialMassFluid.set
    (
        i,
        new dimensionedScalar
        (
            "initialMass",
            dimMass,
            fvc::domainIntegrate
            (
                // thermoFluid[i].alpha2()
                // *thermoFluid[i].thermo2().rho()
                rhoFluid[i]
            ).value()
        )
    );

    Info<< "    creating Phase Change Model \n" << endl;
    PhaseChangeModelFluid.set
    (
        i,
        new singleComponentPhaseChange
        (
            thermoFluid[i],
            thermoFluid[i].p(),
            turbulenceFluid[i],
            advectionSchemesFluid[i].surf()
        )
    );

    Info<< "    Adding to KFluid\n" << endl;
    KFluid.set
    (
        i,
        new volScalarField
        (
            "K",
            0.5*magSqr(UFluid[i])
        )
    );

    fvOptionsFluid.set
    (
        i,
        new fv::options
        (
            fluidRegions[i]
        )
    );

    Info<< "    Adding MRF\n" << endl;
    MRFfluid.set
    (
        i,
        new IOMRFZoneList(fluidRegions[i])
    );

    const dictionary& pimpleDict =
        fluidRegions[i].solutionDict().subDict("PIMPLE");
    pimpleDict.readIfPresent("frozenFlow", frozenFlowFluid[i]);

}
