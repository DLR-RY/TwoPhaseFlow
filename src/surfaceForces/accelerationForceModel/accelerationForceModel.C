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

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(accelerationForceModel, 0);
    defineRunTimeSelectionTable(accelerationForceModel, components);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::accelerationForceModel::accelerationForceModel
(
    const word& type,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    accf_
    (
        IOobject
        (
            "ghf",
            mesh.time().timeName(),
            mesh
        ),
        mesh,
        dimensionedScalar("accf_", dimAcceleration*dimLength, 0.0)
    ),
    acc_
    (
        IOobject
        (
            "gh",
            mesh.time().timeName(),
            mesh
        ),
        mesh,
        dimensionedScalar("acc_", dimAcceleration*dimLength, 0.0)
    )
{

}


// * * * * * * * * * * * * * * Public Access Member Functions  * * * * * * * //


void Foam::accelerationForceModel::calculateAcc()
{
    notImplemented("bool Foam::accelerationForceModel::calculateAcc()");;
}



// ************************************************************************* //
