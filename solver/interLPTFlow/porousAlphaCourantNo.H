/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2021 Johan Roenby
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

Global
    porousAlphaCourantNo

Description
    Calculates and outputs the mean and maximum Courant Numbers.

\*---------------------------------------------------------------------------*/

scalar maxAlphaCo = runTime.controlDict().get<scalar>("maxAlphaCo");
scalar alphaCoNum = 0;
scalar meanAlphaCoNum = 0;

if (mesh.nInternalFaces())
{
    scalarField sumPhi
    (
        mixture.nearInterface()().primitiveField()
       *fvc::surfaceSum(mag(phi))().primitiveField()
    );

    if (porosityEnabled)
    {
        const volScalarField& porosity = tporosity.cref();

        alphaCoNum =
            0.5*gMax(sumPhi/mesh.V().field()/porosity)*runTime.deltaTValue();

        meanAlphaCoNum =
            0.5*(gSum(sumPhi)/gSum(mesh.V().field()*porosity))
           *runTime.deltaTValue();
    }
    else
    {
        alphaCoNum = 0.5*gMax(sumPhi/mesh.V().field())*runTime.deltaTValue();

        meanAlphaCoNum =
            0.5*(gSum(sumPhi)/gSum(mesh.V().field()))*runTime.deltaTValue();
    }
}

Info<< "Interface Courant Number mean: " << meanAlphaCoNum
    << " max: " << alphaCoNum << endl;

// ************************************************************************* //
