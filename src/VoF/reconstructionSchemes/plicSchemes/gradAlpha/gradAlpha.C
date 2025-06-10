/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2019 OpenCFD Ltd.
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

#include "gradAlpha.H"
#include "fvc.H"
#include "leastSquareGrad.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace reconstruction
{
    defineTypeNameAndDebug(gradAlpha, 0);
    addToRunTimeSelectionTable(reconstructionSchemes, gradAlpha, components);
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::reconstruction::gradAlpha::gradSurf(const volScalarField& phi)
{
    addProfilingInFunction(geometricVoF);
    leastSquareGrad<scalar> lsGrad("polyDegree1",mesh_.geometricD());

    zoneDistribute& exchangeFields = zoneDistribute::New(mesh_);

    exchangeFields.setUpCommforZone(interfaceCell_,true);

    Map<vector> mapCC
    (
        exchangeFields.getDatafromOtherProc(interfaceCell_, mesh_.C())
    );
    Map<scalar> mapPhi
    (
        exchangeFields.getDatafromOtherProc(interfaceCell_, phi)
    );

    DynamicField<vector> cellCentre(100);
    DynamicField<scalar> phiValues(100);

    const labelListList& stencil = exchangeFields.getStencil();

    forAll(interfaceLabels_, i)
    {
        const label celli = interfaceLabels_[i];

        cellCentre.clear();
        phiValues.clear();

        for (const label gblIdx : stencil[celli])
        {
            cellCentre.append
            (
                exchangeFields.getValue(mesh_.C(), mapCC, gblIdx)
            );
            phiValues.append
            (
                exchangeFields.getValue(phi, mapPhi, gblIdx)
            );
        }

        cellCentre -= mesh_.C()[celli];
        interfaceNormal_[i] = lsGrad.grad(cellCentre, phiValues);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::reconstruction::gradAlpha::gradAlpha
(
    volScalarField& alpha1,
    const surfaceScalarField& phi,
    const volVectorField& U,
    const dictionary& dict
)
:
    reconstructionSchemes
    (
        typeName,
        alpha1,
        phi,
        U,
        dict
    ),
    mesh_(alpha1.mesh()),
    interfaceNormal_(fvc::grad(alpha1)),
    isoFaceTol_(modelDict().lookupOrDefault<scalar>("isoFaceTol", 1e-8)),
    surfCellTol_(modelDict().lookupOrDefault<scalar>("surfCellTol", 1e-8)),
    sIterPLIC_(mesh_,surfCellTol_)
{
    reconstruct();
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::reconstruction::gradAlpha::reconstruct(bool forceUpdate)
{
    addProfilingInFunction(geometricVoF);
    const bool uptodate = alreadyReconstructed(forceUpdate);

    if (uptodate && !forceUpdate)
    {
        return;
    }

    if (mesh_.topoChanging())
    {
        // Introduced resizing to cope with changing meshes
        if (interfaceCell_.size() != mesh_.nCells())
        {
            interfaceCell_.resize(mesh_.nCells());
        }
    }
    interfaceCell_ = false;

    interfaceLabels_.clear();

    forAll(alpha1_, celli)
    {
        if (sIterPLIC_.isASurfaceCell(alpha1_[celli]))
        {
            interfaceCell_[celli] = true; // is set to false earlier
            interfaceLabels_.append(celli);
        }
    }
    interfaceNormal_.resize(interfaceLabels_.size());
    centre_ = dimensionedVector("centre", dimLength, Zero);
    normal_ = dimensionedVector("normal", dimArea, Zero);

    gradSurf(alpha1_);

    forAll(interfaceLabels_, i)
    {
        const label celli = interfaceLabels_[i];
        if (mag(interfaceNormal_[i]) == 0)
        {
            continue;
        }

        sIterPLIC_.vofCutCell
        (
            celli,
            alpha1_[celli],
            isoFaceTol_,
            100,
            interfaceNormal_[i]
        );

        if (sIterPLIC_.cellStatus() == 0)
        {
            normal_[celli] = sIterPLIC_.surfaceArea();
            centre_[celli] = sIterPLIC_.surfaceCentre();
            if (mag(normal_[celli]) == 0)
            {
                normal_[celli] = Zero;
                centre_[celli] = Zero;
            }
        }
        else
        {
            normal_[celli] = Zero;
            centre_[celli] = Zero;
        }
    }
}


void Foam::reconstruction::gradAlpha::mapAlphaField() const
{
    addProfilingInFunction(geometricVoF);
    // Without this line, we seem to get a race condition
    mesh_.C();

    cutCellPLIC cutCell(mesh_);

    forAll(normal_, celli)
    {
        if (mag(normal_[celli]) != 0)
        {
            vector n = normal_[celli]/mag(normal_[celli]);
            scalar cutValue = (centre_[celli] - mesh_.C()[celli]) & (n);
            cutCell.calcSubCell
            (
                celli,
                cutValue,
                n
            );
            alpha1_[celli] = cutCell.VolumeOfFluid();
        }
    }

    alpha1_.correctBoundaryConditions();
    alpha1_.oldTime () = alpha1_;
    alpha1_.oldTime().correctBoundaryConditions();
}


// ************************************************************************* //
