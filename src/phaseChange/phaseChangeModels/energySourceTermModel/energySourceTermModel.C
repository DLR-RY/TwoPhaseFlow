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
#include "zeroGradientFvPatchFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(energySourceTermModel, 0);
    defineRunTimeSelectionTable(energySourceTermModel, components);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::energySourceTermModel::energySourceTermModel
(
    const word& type,
    const phaseModel& phase1,
    const phaseModel& phase2,
    const compressibleInterPhaseTransportModel& turbModel,
    const volScalarField& p,
    singleComponentSatProp& satModel,
    reconstructionSchemes& surf,
    const dictionary& dict
)
:
    dictionary(dict),
    energySourceTermModelCoeffs_(optionalSubDict(type + "Coeffs")),
    phase1_(phase1),
    phase2_(phase2),
    turbModel_(turbModel),
    p_(p),
    satModel_(satModel),
    surf_(surf)
{

}

// * * * * * * * * * * * * * * Public Access Member Functions  * * * * * * * //


const Foam::dictionary&
Foam::energySourceTermModel::modelDict() const
{
    return energySourceTermModelCoeffs_;
}

Foam::dictionary&
Foam::energySourceTermModel::modelDict()
{
    return energySourceTermModelCoeffs_;
}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //



// ************************************************************************* //
