volScalarField::Internal Sp
(
    IOobject
    (
        "Sp",
        runTime.timeName(),
        mesh
    ),
    mesh,
    dimensionedScalar(dgdt.dimensions(), Zero)
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
    dimensionedScalar(dgdt.dimensions(), Zero)
);

forAll(dgdt, celli)
{
    if (dgdt[celli] > 0.0 && alpha1[celli] > 0.0)
    {
        Sp[celli] -= dgdt[celli]*alpha1[celli];
        Su[celli] += dgdt[celli]*alpha1[celli];
    }
    else if (dgdt[celli] < 0.0 && alpha1[celli] < 1.0)
    {
        Sp[celli] += dgdt[celli]*(1.0 - alpha1[celli]);
    }
}

volScalarField::Internal divU
(
    mesh.moving()
  ? fvc::div(phi + mesh.phi())
  : fvc::div(phi)
);
