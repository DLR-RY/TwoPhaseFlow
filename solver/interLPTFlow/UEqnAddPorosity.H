// Including porosity effects in UEqn following:
// Jensen, B., Jacobsen, N. G., & Christensen, E. D. (2014).
// Investigations on the porous media equations and resistance
// coefficients for coastal structures. Coastal Engineering, 84, 56-72.

if (porosityEnabled)
{
    const volScalarField& porosity = tporosity.cref();

    const word porosityModel("JensenEtAl2014");
    const dictionary& dict =
        porosityProperties.subDict(porosityModel + "Coeffs");
    const dimensionedScalar alpha(dimless/dimArea, dict.get<scalar>("alpha"));
    const dimensionedScalar beta(dimless/dimLength, dict.get<scalar>("beta"));
    const dimensionedScalar d50(dimless, dict.get<scalar>("d50"));
    const dimensionedScalar KC(dimless, dict.get<scalar>("KC"));

    // Generating Darcy-Forchheimer coefficient: F = rho*U*(a + b*|U|)
    // Shoud it be mu or muEff in the equation below?
    {
        // Darcy term
        volScalarField DarcyForchheimerCoeff
        (
            alpha*sqr(1 - porosity)*mixture.mu()/sqr(porosity)/sqr(d50)
        );

        // Adding Forchheimer term
        DarcyForchheimerCoeff += rho*mag(U)
            *beta*(1 + pos(KC)*7.5/KC)*(1 - porosity)/sqr(porosity)/d50;

        // Adding Darcy-Forchheimer term as implicit source term
        UEqn += fvm::Sp(DarcyForchheimerCoeff, U);
    }

    {
        // Generating added mass force coefficient
        const dimensionedScalar gamma_p(dimless, dict.get<scalar>("gamma_p"));
        const volScalarField Cm(gamma_p*(1 - porosity));

        UEqn += Cm*fvm::ddt(rho, U);
        UEqn += Cm*MRF.DDt(rho, U);
    }

    // Dividing both matrix entries and source term by porosity to compensate
    // for the fact that the FVM cell volume averages use division by cell
    // volume V whereas only the cell pore volume, porosity*V, is accessible.
    UEqn *= scalar(1)/porosity;
}
