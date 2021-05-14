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
#include "messageStream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::reconstructionSchemes>
Foam::reconstructionSchemes::New
(
    volScalarField& alpha1,
    const surfaceScalarField& phi,
    const volVectorField& U,
    const dictionary& dict
)
{
    word schemeType("isoAlpha");

    if (!dict.readIfPresent("reconstructionScheme", schemeType))
    {
        Warning
            << "Entry '"
            << "reconstructionScheme" << "' not found in dictionary "
            << dict.name() << nl
            << "using default" << nl;
    }

    Info<< "Selecting reconstructionScheme: " << schemeType << endl;

    auto cstrIter = componentsConstructorTablePtr_->cfind(schemeType);

    if (!cstrIter.found())
    {
        FatalErrorInFunction
            << "Unknown reconstructionSchemes type "
            << schemeType << nl << nl
            << "Valid  reconstructionSchemess are : " << nl
            << componentsConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<reconstructionSchemes>(cstrIter()(alpha1, phi, U, dict));
}


// ************************************************************************* //
