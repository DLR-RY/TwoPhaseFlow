
    scalar cutOff = 1e-3;

    surfaceScalarField kappaLEff
    (
        "kappaLEff",
        fvc::interpolate(mixture.thermo1().kappaEff(turbulence.alphat()))
    );
    surfaceScalarField kappaVEff
    (
        "kappaVEff",
        fvc::interpolate(mixture.thermo2().kappaEff(turbulence.alphat()))
    );

    const labelList& owner = mesh.owner();
    const labelList& neighbour = mesh.neighbour();
    const volScalarField& TSat = phaseChange.satProp().TSat();

    forAll(kappaLEff,faceI)
    {
        const label& own = owner[faceI];
        const label& nei = neighbour[faceI];
        const scalar alphaOwn = alpha1[own];
        const scalar alphaNei = alpha1[nei];

        if(alphaOwn < cutOff || alphaNei < cutOff) // liquid phase
        {
            kappaLEff[faceI] = 0.0;
        }
        if(alphaOwn > (1-cutOff) || alphaNei > (1-cutOff)) // gas phase
        {
            kappaVEff[faceI] = 0.0;
        }
    }

    // energy equation
    fvScalarMatrix TLEqn
    (
        fvm::ddt(rhoCpL,TL)
        + fvm::div(rhoCpLPhi,TL)
	    - fvm::Sp(fvc::ddt(rhoCpL) + fvc::div(rhoCpLPhi), TL)
        - fvm::laplacian(kappaLEff,TL)
        + phaseChange.TSource1()
            ==
        fvOptions(rhoCpL, TL)
    );

    // energy equation
    fvScalarMatrix TVEqn
    (
        fvm::ddt(rhoCpV,TV)
        + fvm::div(rhoCpVPhi,TV)
	    - fvm::Sp(fvc::ddt(rhoCpV) + fvc::div(rhoCpVPhi), TV)
        - fvm::laplacian(kappaVEff,TV)
        - dpdt*alpha2
        + phaseChange.TSource2()
            ==
        fvOptions(rhoCpV, TV)
    );

    label liquidCellSize = 0;
    label gasCellSize = 0;

    labelList FluidStatus(mesh.nCells(),0);

    forAll(alpha1,cellI)
    {
        if (alpha1.oldTime()[cellI] >= (1-cutOff))
        {
            FluidStatus[cellI]= 1;
            liquidCellSize++;
        }
        else if (alpha1.oldTime()[cellI] <= cutOff)
        {
            FluidStatus[cellI]= -1;
            gasCellSize++;
        }
        else
        {
            FluidStatus[cellI]= 0;
        }
    }

    labelList  GasCells(gasCellSize);
    scalarList GasSatTemp(gasCellSize);

    labelList  LiquidCells(liquidCellSize);
    scalarList LiquidSatTemp(liquidCellSize);

    label liquidCellCounter = 0;
    label gasCellCounter = 0;

    forAll(FluidStatus,cellI)
    {
        if (FluidStatus[cellI]== -1) // gas
        {
            GasCells[gasCellCounter] = cellI;
            GasSatTemp[gasCellCounter] = TSat[cellI];
            gasCellCounter++;
        }

        if (FluidStatus[cellI] == 1) // liquid
        {
            LiquidCells[liquidCellCounter]= cellI;
            LiquidSatTemp[liquidCellCounter] = TSat[cellI];
            liquidCellCounter++;
        }
    }

    TL = max(min(TL,TMax),TMin);
    TV = max(min(TV,TMax),TMin);

    TLEqn.setValues(GasCells, GasSatTemp); // gas cells are set to sat Temp
    TVEqn.setValues(LiquidCells, LiquidSatTemp); // liquid cells are set to sat temp

    TL.correctBoundaryConditions();
    TV.correctBoundaryConditions();

    TLEqn.relax();
    TVEqn.relax();

    TLEqn.solve();
    TVEqn.solve();

    Info << "min(TL) " << gMin(TL) << "  max(TL) " << gMax(TL)  << endl;
    Info << "min(TV) " << gMin(TV) << "  max(TV) " << gMax(TV) << endl;

    TL = max(min(TL,TMax),TMin);
    TV = max(min(TV,TMax),TMin);

    T = TL*alpha1 + TV*alpha2;

    T.correctBoundaryConditions();

    mixture.correctThermo(TL,TV);

    mixture.correct();


    Info<< "min(T) " << min(T).value() << endl;
