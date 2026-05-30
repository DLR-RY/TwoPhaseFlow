Info<< " Reading electrical properties\n" << nl << endl;

IOdictionary electricalProperties
(
    IOobject
    (
        "electricalProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

dimensionedScalar permittivity1
    (
        "permittivity1",
        dimensionSet(-1, -3, 4, 0, 0, 2, 0),
        electricalProperties
    );
dimensionedScalar permittivity2
    (
        "permittivity2",
        dimensionSet(-1, -3, 4, 0, 0, 2, 0),
        electricalProperties
    );
    dimensionedScalar ionConductivity1
    (
        "ionConductivity1",
        dimensionSet(-1, -3, 3, 0, 0, 2, 0),
        electricalProperties
    );
 dimensionedScalar ionConductivity2
    (
        "ionConductivity2",
        dimensionSet(-1, -3, 3, 0, 0, 2, 0),
        electricalProperties
    );

Info<< "Reading field chargeDensity\n" << endl;

volScalarField chargeDensity
(
    IOobject
    (
        "chargeDensity",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field potentialDifference\n" << endl;
volScalarField potentialDifference
(
    IOobject
    (
        "potentialDifference",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);


Info<< "Reading field ionConductivity\n" << endl;


volScalarField ionConductivity
(
    IOobject
    (
        "ionConductivity",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    (alpha1*ionConductivity1)+(alpha2*ionConductivity2)
);
// ionConductivity.oldTime();

surfaceScalarField ionConductivityF
(
    IOobject
    (
        "ionConductivityF",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    (fvc::interpolate(alpha1)*ionConductivity1 + fvc::interpolate(alpha2)*ionConductivity2)
);

Info<< "Reading field permittivity\n" << endl;
volScalarField permittivity
(
    IOobject
    (
        "permittivity",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    (alpha1*permittivity1)+(alpha2*permittivity2)
);
// permittivity.oldTime();

surfaceScalarField permittivityF
(
    IOobject
    (
        "permittivityF",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    (fvc::interpolate(alpha1)*permittivity1 + fvc::interpolate(alpha2)*permittivity2)
);

volVectorField electricField
(
    IOobject
    (
        "electricField",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedVector("electricField", dimensionSet(1,1,-3,0,0,-1,0), vector::zero)
);

Info<< "Calculating field maxwellTensor\n" << endl;
volTensorField maxwellTensor
(
    IOobject
    (
        "maxwellTensor",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    permittivity*((electricField*electricField)-(double(0.5)*magSqr(electricField)*tensor::I))
);

volVectorField ehdForce
(
    IOobject
    (
        "ehdForce",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedVector("ehdForce", dimensionSet(1,-2,-2,0,0,0,0), vector::zero)
);


const dictionary& ehdDict =
    mesh.solutionDict().subDict("ehdSolution");

Switch coupleConvChargeDensity
(
    ehdDict.lookupOrDefault<Switch>("coupleConvChargeDensity", true)
);

Switch coupleMaxwellTensor
(
    ehdDict.lookupOrDefault<Switch>("coupleMaxwellTensor", true)
);

label ehdCorr
(
    ehdDict.lookupOrDefault<label>("ehdCorr", 1)
);

word ionConductivityAvg
(
    ehdDict.lookupOrDefault<word>("ionConductivityAvg", "arithmetic")
);

word permittivityAvg
(
    ehdDict.lookupOrDefault<word>("permittivityAvg", "arithmetic")
);

#include "casesEHD.H"



