if(oCorr == 0 && closedVolume) // correct the pressure with the total mass
{
    const polyBoundaryMesh& pbm = mesh.boundaryMesh();
    dimensionedScalar intPsi0 = fvc::domainIntegrate(phaseChange.psi0());
    const fvPatchList& patches = mesh.boundary();

    forAll(patches, patchI)
    {
        const polyPatch& pp = pbm[patchI];
        if
        (
            !isType<processorFvPatch>(mesh.boundary()[patchI]) &&
            !isA<emptyPolyPatch>(pp) &&
            !isA<wedgePolyPatch>(pp)
        )
        {
            dimensionedScalar massFluxP
            (
                "0",
                dimensionSet(1, 0, -1, 0, 0, 0, 0),
                gSum(rhoPhi.boundaryField()[patchI])
            );
            // patch vector is always pointing outwards -> minus
            initialMass -= massFluxP*runTime.deltaT();
        }
    }

}
