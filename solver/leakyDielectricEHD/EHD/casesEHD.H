Info << "Coupling Maxwell tensor with UEqn: " << coupleMaxwellTensor << endl;

Info << "EHD correction loop: " << ehdCorr <<endl;

if (ionConductivityAvg == "arithmetic")
{
    Info << "Using arithmetic averaging for ionConductivityAvg" << endl;
}
else if (ionConductivityAvg == "harmonic")
{
    Info << "Using harmonic averaging for ionConductivityAvg" << endl;
}
else
{
    FatalErrorInFunction
        << "Unknown interfaceAvg for ionConductivityAvg: " << ionConductivityAvg << nl
        << "Valid options are: arithmetic, harmonic"
        << exit(FatalError);
}

if (permittivityAvg == "arithmetic")
{
    Info << "Using arithmetic averaging for permittivityAvg" << endl;
}
else if (permittivityAvg == "harmonic")
{
    Info << "Using harmonic averaging for permittivityAvg" << endl;
}
else
{
    FatalErrorInFunction
        << "Unknown interfaceAvg for permittivityAvg: " << permittivityAvg << nl
        << "Valid options are: arithmetic, harmonic"
        << exit(FatalError);
}

