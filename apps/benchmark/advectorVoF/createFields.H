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

Info<< "Reading transportProperties\n" << endl;
immiscibleIncompressibleTwoPhaseMixture mixture(U, phi);

volScalarField& alpha1(mixture.alpha1());
volScalarField& alpha2(mixture.alpha2());

const dimensionedScalar& rho1 = mixture.rho1();
const dimensionedScalar& rho2 = mixture.rho2();

scalar V0 = fvc::domainIntegrate(alpha1).value();

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

IOdictionary fvSolutionDict
(
    IOobject
    (
        "fvSolution",
        mesh.time().system(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);


autoPtr<functionObjects::setFlow> flow;

if (fvSolutionDict.found("setFlow"))
{
    dictionary setFlowDict = fvSolutionDict.subDict("setFlow");
    flow.reset(new functionObjects::setFlow("setFlow",runTime,setFlowDict));
}

const dictionary& isoAdvectorDict = mesh.solutionDict().subDict("isoAdvector");
const scalar period = isoAdvectorDict.lookupOrDefault<scalar>("period", 0.0);
scalar reverseTime = isoAdvectorDict.lookupOrDefault<scalar>("reverseTime", 0.0);
scalar spirallingFlow = isoAdvectorDict.lookupOrDefault<scalar>("spirallingFlow", 0.0);

//Copies of initial U and phi for use when flow is periodic
surfaceScalarField phi0 = phi;
volVectorField U0 = U;
