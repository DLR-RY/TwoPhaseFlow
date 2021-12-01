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

#include "macroModel.H"
#include "zeroGradientFvPatchFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(macroModel, 0);
    defineRunTimeSelectionTable(macroModel, components);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::macroModel::macroModel
(
    const word& type,
    const phaseModel& phase1,
    const phaseModel& phase2,
    const volScalarField& p,
    singleComponentSatProp& satModel,
    const compressibleInterPhaseTransportModel& turbModel,
    const dictionary& dict
)
:
    dictionary(dict),
    macroModelCoeffs_(optionalSubDict(type + "Coeffs")),
    phase1_(phase1),
    phase2_(phase2),
    p_(p),
    satModel_(satModel),
    turbModel_(turbModel)
{

}

// * * * * * * * * * * * * * * Public Access Member Functions  * * * * * * * //


const Foam::dictionary&
Foam::macroModel::modelDict() const
{
    return macroModelCoeffs_;
}

Foam::dictionary&
Foam::macroModel::modelDict()
{
    return macroModelCoeffs_;
}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //



// ************************************************************************* //
