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
    porousCourantNo

Description
    Calculates and outputs the mean and maximum Courant Numbers.

\*---------------------------------------------------------------------------*/

scalar CoNum = 0;
scalar meanCoNum = 0;

{
    scalarField sumPhi
    (
        fvc::surfaceSum(mag(phi))().primitiveField()
    );

    if (porosityEnabled)
    {
        const volScalarField& porosity = tporosity.cref();

        CoNum =
            0.5*gMax(sumPhi/mesh.V().field()/porosity)
           *runTime.deltaTValue();

        meanCoNum =
            0.5*(gSum(sumPhi)/gSum(mesh.V().field()*porosity))
           *runTime.deltaTValue();
    }
    else
    {
        CoNum = 0.5*gMax(sumPhi/mesh.V().field())*runTime.deltaTValue();

        meanCoNum =
            0.5*(gSum(sumPhi)/gSum(mesh.V().field()))*runTime.deltaTValue();
    }
}

Info<< "Courant Number mean: " << meanCoNum
    << " max: " << CoNum << endl;

// ************************************************************************* //
