/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019-2020 DLR
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

#include "reconstructionSchemes.H"
#include "cutCellPLIC.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(reconstructionSchemes, 0);
    defineRunTimeSelectionTable(reconstructionSchemes, components);
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

bool Foam::reconstructionSchemes::alreadyReconstructed(bool forceUpdate) const
{
    const Time& runTime = alpha1_.mesh().time();

    label& curTimeIndex = timeIndexAndIter_.first();
    label& curIter = timeIndexAndIter_.second();

    // Reset timeIndex and curIter
    if (curTimeIndex < runTime.timeIndex())
    {
        curTimeIndex = runTime.timeIndex();
        curIter = 0;
        return false;
    }

    if (forceUpdate)
    {
        curIter = 0;
        return false;
    }

    // Always reconstruct when subcycling
    if (runTime.subCycling() != 0)
    {
        return false;
    }

    ++curIter;
    if (curIter > 1)
    {
        return true;
    }

    return false;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::reconstructionSchemes::reconstructionSchemes
(
    const word& type,
    volScalarField& alpha1,
    const surfaceScalarField& phi,
    const volVectorField& U,
    const dictionary& dict
)
:
    IOdictionary
    (
        IOobject
        (
            "reconstructionScheme",
            alpha1.time().constant(),
            alpha1.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        )
    ),
    reconstructionSchemesCoeffs_(dict),
    alpha1_(alpha1),
    phi_(phi),
    U_(U),
    normal_
    (
        IOobject
        (
            IOobject::groupName("interfaceNormal", alpha1.group()),
            alpha1_.mesh().time().timeName(),
            alpha1_.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        alpha1_.mesh(),
        dimensionedVector(dimArea, Zero)
    ),
    centre_
    (
        IOobject
        (
            IOobject::groupName("interfaceCentre", alpha1.group()),
            alpha1_.mesh().time().timeName(),
            alpha1_.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        alpha1_.mesh(),
        dimensionedVector(dimLength, Zero)
    ),
    interfaceCell_(alpha1_.mesh().nCells(), false),
    interfaceLabels_(0.2*alpha1_.mesh().nCells()),
    timeIndexAndIter_(0, 0)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //



// ************************************************************************* //
