/*********** Solving Herrera Model ************/

// Inspired from Basilisk method of equation solving
fvScalarMatrix potentialDifferenceEqn
(
    fvm::laplacian(ionConductivityF, potentialDifference)
);
// potentialDifferenceEqn.relax();
potentialDifferenceEqn.solve();

// Calculating electric Field
electricField = -fvc::grad(potentialDifference);

//(LopezHerrera2011, Roghair2015)
chargeDensity = -fvc::laplacian(permittivityF, potentialDifference);

permittivity = fvc::average(permittivityF);
ionConductivity = fvc::average(ionConductivityF);

maxwellTensor = permittivity*((electricField*electricField)-(double(0.5)*magSqr(electricField)*tensor::I));

if (coupleMaxwellTensor)
{
    ehdForce =  fvc::div(maxwellTensor);
}
