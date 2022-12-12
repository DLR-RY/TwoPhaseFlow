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

#include "isoAlpha.H"
#include "addToRunTimeSelectionTable.H"
#include "cutCellPLIC.H"
#include "profiling.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace reconstruction
{
    defineTypeNameAndDebug(isoAlpha, 0);
    addToRunTimeSelectionTable(reconstructionSchemes,isoAlpha, components);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::reconstruction::isoAlpha::isoAlpha
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
    // Interpolation data
    ap_(mesh_.nPoints()),

    // Tolerances and solution controls
    isoFaceTol_(modelDict().lookupOrDefault<scalar>("isoFaceTol", 1e-8)),
    surfCellTol_(modelDict().lookupOrDefault<scalar>("surfCellTol", 1e-8)),
    sIterIso_(mesh_,ap_,surfCellTol_)
{
    reconstruct();
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::reconstruction::isoAlpha::reconstruct(bool forceUpdate)
{
    addProfilingInFunction(geometricVoF);
    const bool uptodate = alreadyReconstructed(forceUpdate);

    if (uptodate && !forceUpdate)
    {
        return;
    }

    // Interpolating alpha1 cell centre values to mesh points (vertices)
    if (mesh_.topoChanging())
    {
        // Introduced resizing to cope with changing meshes
        if (ap_.size() != mesh_.nPoints())
        {
            ap_.resize(mesh_.nPoints());

        }
        if (interfaceCell_.size() != mesh_.nCells())
        {
            interfaceCell_.resize(mesh_.nCells());
        }
    }
    ap_ = volPointInterpolation::New(mesh_).interpolate(alpha1_);

    DynamicList<List<point>> facePts;

    interfaceLabels_.clear();

    forAll(alpha1_,cellI)
    {
        if (sIterIso_.isASurfaceCell(alpha1_[cellI]))
        {
            interfaceLabels_.append(cellI);

            sIterIso_.vofCutCell
            (
                cellI,
                alpha1_[cellI],
                isoFaceTol_,
                100
            );

            if (sIterIso_.cellStatus() == 0)
            {
                normal_[cellI] = sIterIso_.surfaceArea();
                centre_[cellI] = sIterIso_.surfaceCentre();
                if (mag(normal_[cellI]) != 0)
                {
                    interfaceCell_[cellI] = true;
                }
                else
                {
                    interfaceCell_[cellI] = false;
                }
            }
            else
            {
                normal_[cellI] = Zero;
                centre_[cellI] = Zero;
                interfaceCell_[cellI] = false;
            }
         }
         else
         {
            normal_[cellI] = Zero;
            centre_[cellI] = Zero;
            interfaceCell_[cellI] = false;
         }
    }
}


// ************************************************************************* //
