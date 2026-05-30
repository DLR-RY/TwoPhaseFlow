volScalarField::Internal Sp
(
    IOobject
    (
        "Sp",
        runTime.timeName(),
        mesh
    ),
    mesh,
    dimensionedScalar(dimless/dimTime, Zero)
);

volScalarField::Internal Su
(
    IOobject
    (
        "Su",
        runTime.timeName(),
        mesh
    ),
    mesh,
    dimensionedScalar(dimless/dimTime, Zero)
);

volScalarField::Internal divU
(
    mesh.moving()
  ? fvc::div(phi + mesh.phi())
  : fvc::div(phi)
);
