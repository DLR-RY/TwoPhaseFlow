/*---------------------------------------------------------------------------*\
            Copyright (c) 2017-2019, German Aerospace Center (DLR)
-------------------------------------------------------------------------------
License
    This file is part of the VoFLibrary source code library, which is an
	unofficial extension to OpenFOAM.
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

#include "energySourceTermModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::energySourceTermModel>
Foam::energySourceTermModel::New
(
    const phaseModel& phase1,
    const phaseModel& phase2,
    const compressibleInterPhaseTransportModel& turbModel,
    const volScalarField& p,
    singleComponentSatProp& satModel,
    reconstructionSchemes& surf,
    const dictionary& dict
)
{
    word energySourceTermModelTypeName
    (
        dict.lookup("energySourceTermModel")
    );

    Info<< "Selecting energySourceTermModel model "
        << energySourceTermModelTypeName << endl;


    auto ctorPtr = componentsConstructorTable(energySourceTermModelTypeName);

    if (!ctorPtr)
    {
        FatalErrorInLookup
        (
            "energySourceTermModel",
            energySourceTermModelTypeName,
            *componentsConstructorTablePtr_
        ) << exit(FatalError);
    }

    return autoPtr<energySourceTermModel>
    (
        ctorPtr
        (
            phase1,
            phase2,
            turbModel,
            p,
            satModel,
            surf,
            dict
        )
    );
}


// ************************************************************************* //
