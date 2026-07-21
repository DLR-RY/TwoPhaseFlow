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

#include "isoSurface.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace reconstruction
{
    defineTypeNameAndDebug(isoSurface, 0);
    addToRunTimeSelectionTable(reconstructionSchemes,isoSurface, components);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::reconstruction::isoSurface::isoSurface
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
    cutCell_(mesh_,ap_),

    // Tolerances and solution controls
    iso_(modelDict().lookupOrDefault<scalar>("iso", 0.5))
{
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::reconstruction::isoSurface::~isoSurface()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::reconstruction::isoSurface::reconstruct(bool forceUpdate)
{
    const bool uptodate = alreadyReconstructed(forceUpdate);

    if (uptodate && !forceUpdate)
    {
        return;
    }
    // Interpolating alpha1 cell centre values to mesh points (vertices)
    ap_ = volPointInterpolation::New(mesh_).interpolate(alpha1_);

    forAll(alpha1_,cellI)
    {
            cutCell_.calcSubCell(cellI, iso_);

            if(cutCell_.cellStatus() == 0)
            {
                normal_[cellI] = cutCell_.faceArea();
                centre_[cellI] = cutCell_.faceCentre();
                if(mag(normal_[cellI]) != 0)
                {
                    interfaceCell_[cellI]=true;
                }
                else
                {
                    interfaceCell_[cellI]=false;
                    normal_[cellI] = vector::zero;
                    centre_[cellI] = vector::zero;
                }
            }
            else
            {
                normal_[cellI] = vector::zero;
                centre_[cellI] = vector::zero;
                interfaceCell_[cellI]=false;
            }

    }

}

// ************************************************************************* //
