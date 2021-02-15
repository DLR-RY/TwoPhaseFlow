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

//#include "physicoChemicaClausiusClapeyrons.H"
#include "singleComponentFunction.H"
#include "addToRunTimeSelectionTable.H"
#include "zeroGradientFvPatchFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(singleComponentFunction, 0);
    addToRunTimeSelectionTable(singleComponentSatProp,singleComponentFunction, components);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::singleComponentFunction::singleComponentFunction
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    singleComponentSatProp
    (
        mesh,
        dict
    ),
    pSatFunc_(Function1<scalar>::New("pSat", dict)),
    TSatFunc_(Function1<scalar>::New("TSat", dict)),
    LFunc_(Function1<scalar>::New("L", dict))
{
}

// * * * * * * * * * * * * * * Public Access Member Functions  * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //


Foam::scalar
Foam::singleComponentFunction::satT(scalar p) const
{
    return TSatFunc_->value(p);
}

Foam::scalar
Foam::singleComponentFunction::satP(Foam::scalar T) const
{
    return pSatFunc_->value(T);
}

Foam::scalar
Foam::singleComponentFunction::satL(Foam::scalar p) const
{
    return LFunc_->value(p);
}

// ************************************************************************* //
