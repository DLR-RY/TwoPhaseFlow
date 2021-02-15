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

#include "sampledInterface.H"
#include "volFieldsFwd.H"
#include "pointFields.H"
#include "volPointInterpolation.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::sampledInterface::sampleOnFaces
(
    const interpolation<Type>& sampler
) const
{
    updateGeometry();  // Recreate geometry if time has changed

    return sampledSurface::sampleOnFaces
    (
        sampler,
        surface().meshCells(),
        surface(),
        points()
    );
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::sampledInterface::sampleOnPoints
(
    const interpolation<Type>& interpolator
) const
{
    notImplemented("interpolation on the points values is currently not implemented");
    updateGeometry();  // Recreate geometry if time has changed

    // Assume volPointInterpolation for the point field!
    const auto& volFld = interpolator.psi();

    if (subMeshPtr_.valid())
    {
        auto tvolSubFld = subMeshPtr_().interpolate(volFld);
        const auto& volSubFld = tvolSubFld();

        auto tpointFld =
            volPointInterpolation::New(volSubFld.mesh()).interpolate(volSubFld);

        return surface().interpolate
        (
            volSubFld,
            tpointFld()
        );
    }


    auto tpointFld =
        volPointInterpolation::New(volFld.mesh()).interpolate(volFld);

    return surface().interpolate
    (
        volFld,
        tpointFld()
    );
}


// ************************************************************************* //
