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


#include "hardtWondra.H"
#include "addToRunTimeSelectionTable.H"
#include "zeroGradientFvPatchFields.H"

#include "mathematicalConstants.H"
#include "surfaceInterpolate.H"
#include "fvcDiv.H"
#include "fvcGrad.H"
#include "fvcSnGrad.H"

#include "fvCFD.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(hardtWondra, 0);
    addToRunTimeSelectionTable(massSourceTermModel,hardtWondra, components);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::hardtWondra::hardtWondra
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
    ),
    cutoff_(modelDict().lookupOrDefault<scalar>("cutoff",1e-3)),
    spread_(modelDict().lookupOrDefault<scalar>("spread",3))
{

}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //


Foam::tmp<Foam::volScalarField>
Foam::hardtWondra::massSource( volScalarField& rhoSource)
{
    tmp<volScalarField> massSource(rhoSource * 0.0);
    volScalarField& massSourceRef = massSource.ref();
    const fvMesh& mesh = phase1_.mesh();

    dimensionedScalar DPsi
    (
        "DPsi",
        dimensionSet(0,2,0,0,0,0,0),
        spread_/sqr(gAverage(mesh.nonOrthDeltaCoeffs()))
    );

    dimensionedScalar intPsi0 = fvc::domainIntegrate(rhoSource);

    volScalarField psi
    (
        IOobject
        (
            "psi",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("0", dimDensity/dimTime, 0),
        "zeroGradient"
    );

    //- Smearing of source term field
    fvScalarMatrix psiEqn
    (
        fvm::Sp(scalar(1),psi) - fvm::laplacian(DPsi,psi) == rhoSource
    );

    psiEqn.solve();

    // Cut cells with cutoff < alpha1 < 1-cutoff
    // and rescale remaining source term field
    dimensionedScalar intPsiLiquid
    (
        "intPsiLiquid",
        dimensionSet(1,0,-1,0,0,0,0),
        0.0
    );
    dimensionedScalar intPsiVapor
    (
        "intPsiLiquid",
        dimensionSet(1,0,-1,0,0,0,0),
        0.0
    );

    forAll(mesh.C(),celli)
    {
        if (phase1_[celli] < cutoff_)
        {
            intPsiVapor.value() +=
                (1.0-phase1_[celli])*psi[celli]*mesh.V()[celli];
        }
        else if (phase1_[celli] > 1.0-cutoff_)
        {
            intPsiLiquid.value() +=
                phase1_[celli]*psi[celli]*mesh.V()[celli];
        }
    }

    //- Calculate Nl and Nv
    dimensionedScalar Nl ("Nl", dimless, 2.0);
    dimensionedScalar Nv ("Nv", dimless, 2.0);

    reduce(intPsiLiquid.value(),sumOp<scalar>());
    reduce(intPsiVapor.value(),sumOp<scalar>());

    if (intPsiLiquid.value() > 1e-99)
    {
        Nl = intPsi0/intPsiLiquid;
    }
    if (intPsiVapor.value() > 1e-99)
    {
        Nv = intPsi0/intPsiVapor;
    }


    //- Set source terms in cells with alpha1 < cutoff or alpha1 > 1-cutoff
    forAll(mesh.C(),celli)
    {
        if (phase1_[celli] < cutoff_)
        {
            massSourceRef[celli] = Nv.value()*(1.0-phase1_[celli])*psi[celli];
        }
        else if (phase1_[celli] > 1.0-cutoff_)
        {
            massSourceRef[celli] = -Nl.value()*phase1_[celli]*psi[celli];
        }
        else
        {
            massSourceRef[celli] = 0.0;
        }
    }

   return massSource;

}


Foam::tmp<Foam::volScalarField>
Foam::hardtWondra::alphaSource( volScalarField& rhoSource)
{
    tmp<volScalarField> alphaSource
    (
        rhoSource / phase1_.thermo().rho()
    );

   return alphaSource;
}
