/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "multiCourantNo.H"
#include "fvc.H"

Foam::Pair<Foam::scalar> Foam::multiCourantNo
(
    const fvMesh& mesh,
    const Time& runTime,
    const surfaceScalarField& phi,
    const twoPhaseModelThermo& mixture
)
{
    scalarField sumPhi
    (
        fvc::surfaceSum(mag(phi))().internalField()
    );

    scalar CoNum = 0.5*gMax(sumPhi/mesh.V().field())*runTime.deltaTValue();

    scalar meanCoNum =
        0.5*(gSum(sumPhi)/gSum(mesh.V().field()))*runTime.deltaTValue();


    scalar alphaCoNum = 0.0;
    scalar meanAlphaCoNum = 0.0;

    scalarField sumPhiAlpha
    (
        mixture.nearInterface()().primitiveField()
       *fvc::surfaceSum(mag(phi))().primitiveField()
    );

    alphaCoNum = 0.5*gMax(sumPhiAlpha/mesh.V().field())*runTime.deltaTValue();

    meanAlphaCoNum =
        0.5*(gSum(sumPhiAlpha)/gSum(mesh.V().field()))*runTime.deltaTValue();


    Info<< "Interface Courant Number mean: " << meanAlphaCoNum
        << " max: " << alphaCoNum << endl;

    Info<< "Region: " << mesh.name() << " Courant Number mean: " << meanCoNum
        << " max: " << CoNum << endl;

    Pair<scalar> CoNums(CoNum,alphaCoNum);

    return CoNums;
}

Foam::Field<Foam::scalar> Foam::multiCourantNo
(
    const fvMesh& mesh,
    const Time& runTime,
    const surfaceScalarField& phi,
    const twoPhaseModelThermo& mixture,
    surfaceForces& surfForces
)
{
    const volScalarField& rho1 = mixture.thermo1().rho();
    const volScalarField& rho2 = mixture.thermo2().rho();
    scalarField sumPhi
    (
    fvc::surfaceSum(mag(phi))().internalField()
    );

    scalar CoNum = 0.5*gMax(sumPhi/mesh.V().field())*runTime.deltaTValue();

    scalar meanCoNum =
        0.5*(gSum(sumPhi)/gSum(mesh.V().field()))*runTime.deltaTValue();


    scalar alphaCoNum = 0.0;
    scalar meanAlphaCoNum = 0.0;

    scalarField sumPhiAlpha
    (
        mixture.nearInterface()().primitiveField()
       *fvc::surfaceSum(mag(phi))().primitiveField()
    );

    alphaCoNum = 0.5*gMax(sumPhiAlpha/mesh.V().field())*runTime.deltaTValue();

    meanAlphaCoNum =
        0.5*(gSum(sumPhiAlpha)/gSum(mesh.V().field()))*runTime.deltaTValue();

    scalar capillaryNum = runTime.deltaTValue()/surfForces.capillaryDt(rho1,rho2);


    Info<< "Interface Courant Number mean: " << meanAlphaCoNum
        << " max: " << alphaCoNum << endl;

    Info<< "Region: " << mesh.name() << " Courant Number mean: " << meanCoNum
        << " max: " << CoNum << endl;

    Info<< "Capillary Number: " << capillaryNum << endl;

    Field<scalar> CoNums(3);
    CoNums[0] = CoNum;
    CoNums[1] = alphaCoNum;
    CoNums[2] = capillaryNum;

    return CoNums;
}

// ************************************************************************* //
