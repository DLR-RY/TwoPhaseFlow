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


#include "noPhaseChange.H"
#include "addToRunTimeSelectionTable.H"
#include "zeroGradientFvPatchFields.H"

#include "mathematicalConstants.H"
#include "surfaceInterpolate.H"
#include "fvcDiv.H"
#include "fvcGrad.H"
#include "fvcSnGrad.H"

#include "centredCPCCellToCellStencilObject.H"
#include "SortableList.H"

#include "fvCFD.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(noPhaseChange, 0);
    addToRunTimeSelectionTable(massSourceTermModel,noPhaseChange, components);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::noPhaseChange::noPhaseChange
(
    const phaseModel& phase1,
    const phaseModel& phase2,
    const volScalarField& p,
    singleComponentSatProp& satModel,
    reconstructionSchemes& surf,
    const dictionary& dict
)
:
    massSourceTermModel
    (
        typeName,
        phase1,
        phase2,
        p,
        satModel,
        surf,
        dict
    )

{

}
// * * * * * * * * * * * * * * Protected Access Member Functions  * * * * * * //

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //


// ************************************************************************* //
Foam::tmp<Foam::volScalarField>
Foam::noPhaseChange::massSource( volScalarField& rhoSource)
{
    tmp<volScalarField> massSource(rhoSource * 0.0);
    rhoSource *= 0.0;

    return massSource;
}

Foam::tmp<Foam::volScalarField>
Foam::noPhaseChange::alphaSource( volScalarField& rhoSource)
{
    dimensionedScalar scale("scale", dimless/dimDensity,0.0);
    tmp<volScalarField> alphaSource(rhoSource * scale);

   return alphaSource;
}
