Info<< "\nConstructing coal cloud" << endl;
basicKinematicCloud bomb
(
    "kinematicCloud1",
    rho,
    U,
    mu,
    Foam::vector(0,-9.81,0)
);
