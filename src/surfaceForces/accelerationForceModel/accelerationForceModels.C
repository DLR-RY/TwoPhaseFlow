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

#include "accelerationForceModel.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::accelerationForceModel>
Foam::accelerationForceModel::New
(
    const dictionary& dict,
    const fvMesh& mesh
)
{

    word accelerationForceModelTypeName
    (
        dict.getCompat<word>
        (
            "accelerationForceModel",
            {{"accelerationModel", 1.1}}
        )
    );

    Info<< "Selecting surfaceTension model "
        << accelerationForceModelTypeName << endl;

    auto ctorPtr = componentsConstructorTable(accelerationForceModelTypeName);

    if (!ctorPtr)
    {
        FatalErrorInLookup
        (
            "accelerationForceModel",
            accelerationForceModelTypeName,
            *componentsConstructorTablePtr_
        ) << exit(FatalError);
    }

    return autoPtr<accelerationForceModel>
    (
        ctorPtr
        (
            dict,
            mesh
        )
    );
}


// ************************************************************************* //
