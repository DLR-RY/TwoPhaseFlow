    scalar CoNum = -GREAT;
    scalar alphaCoNum = -GREAT;
    scalar capillaryNum = -GREAT;

    forAll(fluidRegions, regionI)
    {
        Field<scalar> CourantNos =  multiCourantNo
        (
            fluidRegions[regionI],
            runTime,
            phiFluid[regionI],
            thermoFluid[regionI],
            surfForcesFluid[regionI]
        );
        CoNum = max
        (
            CourantNos[0],
            CoNum
        );

        alphaCoNum = max
        (
            CourantNos[1],
            alphaCoNum
        );

        capillaryNum = max
        (
            CourantNos[2],
            capillaryNum
        );
    }
