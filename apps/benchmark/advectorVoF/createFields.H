    Info<< "Reading field U\n" << endl;
    volVectorField U
    (
        IOobject
        (
            "U",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    #include "createPhi.H"

/*    Info<< "Reading field alpha1\n" << endl;
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
    );*/
Info<< "Reading transportProperties\n" << endl;
immiscibleIncompressibleTwoPhaseMixture mixture(U, phi);

volScalarField& alpha1(mixture.alpha1());
volScalarField& alpha2(mixture.alpha2());

    Info<< "Reading field alpha1\n" << endl;
   /* volScalarField alphaIso
    (
        IOobject
        (
            "alphaIso",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        alpha1
    );*/

const dimensionedScalar& rho1 = mixture.rho1();
const dimensionedScalar& rho2 = mixture.rho2();

    //Initial Volume of Fluid (VOF)
    scalar V0 = fvc::domainIntegrate(alpha1).value();
//gSum(mesh.V()*alpha1.internalField()());


    Info << "sum(alpha*V) = " << V0
     << ", 1-max(alpha1) = " << 1 - gMax(alpha1.internalField())
     << "\t min(alpha1) = " << gMin(alpha1.internalField()) << endl;

    //IsoAdvector related stuff
    //isoAdvection advectorIso(alphaIso,phi,U);
    // Construct incompressible turbulence model
    autoPtr<advectionSchemes> advector
    (
        advectionSchemes::New(alpha1,phi,U)
    );

    const dictionary& isoAdvectorDict = mesh.solutionDict().subDict("isoAdvector");
    const scalar period = isoAdvectorDict.lookupOrDefault<scalar>("period", 0.0);
    scalar reverseTime = isoAdvectorDict.lookupOrDefault<scalar>("reverseTime", 0.0);
    scalar spirallingFlow = isoAdvectorDict.lookupOrDefault<scalar>("spirallingFlow", 0.0);

    //Copies of initial U and phi for use when flow is periodic
    surfaceScalarField phi0 = phi;
    volVectorField U0 = U;

    //Optional writing out of surface cell set and bounded cell sets
    IOobject::writeOption writeOpt = IOobject::NO_WRITE;
/*    const bool printSurfCells
        = isoAdvectorDict.lookupOrDefault<bool>("printSurfCells", false);
    if (printSurfCells)
    {
        writeOpt = IOobject::AUTO_WRITE;
    }*/

/*    cellSet surfCells
    (
        IOobject
        (
            "surfCells",
            runTime.timeName(),
            runTime,
            IOobject::NO_READ,
            writeOpt
        )
    );

    writeOpt = IOobject::NO_WRITE;
    const bool printBoundCells
        = isoAdvectorDict.lookupOrDefault<bool>("printBoundCells", false);
    if (printBoundCells)
    {
        writeOpt = IOobject::AUTO_WRITE;
    }
    cellSet boundCells
    (
        IOobject
        (
            "boundCells",
            runTime.timeName(),
            runTime,
            IOobject::NO_READ,
            writeOpt
        )
    );*/
