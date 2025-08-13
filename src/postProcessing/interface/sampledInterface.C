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
#include "dictionary.H"
#include "volFields.H"
#include "volPointInterpolation.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(sampledInterface, 0);
    addNamedToRunTimeSelectionTable
    (
        sampledSurface,
        sampledInterface,
        word,
        interface
    );
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //




bool Foam::sampledInterface::updateGeometry() const
{
    const fvMesh& fvm = static_cast<const fvMesh&>(mesh());



    // No update needed
    if (fvm.time().timeIndex() == prevTimeIndex_)
    {
        return false;
    }

    // Get any subMesh
    if (zoneID_.index() != -1 && !subMeshPtr_.valid())
    {
        const polyBoundaryMesh& patches = mesh().boundaryMesh();

        // Patch to put exposed internal faces into
        const label exposedPatchi = patches.findPatchID(exposedPatchName_);

        if (debug)
        {
            Info<< "Allocating subset of size "
                << mesh().cellZones()[zoneID_.index()].size()
                << " with exposed faces into patch "
                << patches[exposedPatchi].name() << endl;
        }

        subMeshPtr_.reset
        (
            new fvMeshSubset(fvm)
        );
        subMeshPtr_().setCellSubset
        (
            labelHashSet(mesh().cellZones()[zoneID_.index()]),
            exposedPatchi
        );
    }


    prevTimeIndex_ = fvm.time().timeIndex();


    // Clear any stored topo
    surfPtr_.clear();

    // Clear derived data
    clearGeom();

//    if (subMeshPtr_.valid())
//    {
//        const volScalarField& vfld = *volSubFieldPtr_;

//        surfPtr_.reset
//        (
//            new interface
//            (
//                vfld.mesh()
//            )
//        );
//    }
//    else
    {
//        const volScalarField& vfld = *volFieldPtr_;


        surfPtr_.reset
        (
            new interface
            (
                fvm
            )
        );
    }



    return true;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::sampledInterface::sampledInterface
(
    const word& name,
    const polyMesh& mesh,
    const dictionary& dict
)
:
    sampledSurface(name, mesh, dict),
    zoneID_(dict.lookupOrDefault("zone", word::null), mesh.cellZones()),
    exposedPatchName_(word::null),
    surfPtr_(nullptr),
    prevTimeIndex_(-1),
    storedVolFieldPtr_(nullptr),
    volFieldPtr_(nullptr),
    pointFieldPtr_(nullptr)
{

    if (zoneID_.index() != -1)
    {
        dict.lookup("exposedPatchName") >> exposedPatchName_;

        if (mesh.boundaryMesh().findPatchID(exposedPatchName_) == -1)
        {
            FatalIOErrorInFunction
            (
                dict
            )   << "Cannot find patch " << exposedPatchName_
                << " in which to put exposed faces." << endl
                << "Valid patches are " << mesh.boundaryMesh().names()
                << exit(FatalIOError);
        }

        if (debug && zoneID_.index() != -1)
        {
            Info<< "Restricting to cellZone " << zoneID_.name()
                << " with exposed internal faces into patch "
                << exposedPatchName_ << endl;
        }
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::sampledInterface::~sampledInterface()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::sampledInterface::needsUpdate() const
{
    const fvMesh& fvm = static_cast<const fvMesh&>(mesh());

    return fvm.time().timeIndex() != prevTimeIndex_;
}


bool Foam::sampledInterface::expire()
{
    surfPtr_.clear();
    subMeshPtr_.clear();

    // Clear derived data
    clearGeom();

    // already marked as expired
    if (prevTimeIndex_ == -1)
    {
        return false;
    }

    // force update
    prevTimeIndex_ = -1;
    return true;
}


bool Foam::sampledInterface::update()
{
    return updateGeometry();
}

Foam::tmp<Foam::scalarField> Foam::sampledInterface::sample
(
    const interpolation<scalar>& sampler
) const
{
    return sampleOnFaces(sampler);
}


Foam::tmp<Foam::vectorField> Foam::sampledInterface::sample
(
    const interpolation<vector>& sampler
) const
{
    return sampleOnFaces(sampler);
}


Foam::tmp<Foam::sphericalTensorField> Foam::sampledInterface::sample
(
    const interpolation<sphericalTensor>& sampler
) const
{
    return sampleOnFaces(sampler);
}


Foam::tmp<Foam::symmTensorField> Foam::sampledInterface::sample
(
    const interpolation<symmTensor>& sampler
) const
{
    return sampleOnFaces(sampler);
}


Foam::tmp<Foam::tensorField> Foam::sampledInterface::sample
(
    const interpolation<tensor>& sampler
) const
{
    return sampleOnFaces(sampler);
}


Foam::tmp<Foam::scalarField> Foam::sampledInterface::interpolate
(
    const interpolation<scalar>& interpolator
) const
{
    return sampleOnPoints(interpolator);
}


Foam::tmp<Foam::vectorField> Foam::sampledInterface::interpolate
(
    const interpolation<vector>& interpolator
) const
{
    return sampleOnPoints(interpolator);
}

Foam::tmp<Foam::sphericalTensorField> Foam::sampledInterface::interpolate
(
    const interpolation<sphericalTensor>& interpolator
) const
{
    return sampleOnPoints(interpolator);
}


Foam::tmp<Foam::symmTensorField> Foam::sampledInterface::interpolate
(
    const interpolation<symmTensor>& interpolator
) const
{
    return sampleOnPoints(interpolator);
}


Foam::tmp<Foam::tensorField> Foam::sampledInterface::interpolate
(
    const interpolation<tensor>& interpolator
) const
{
    return sampleOnPoints(interpolator);
}



void Foam::sampledInterface::print(Ostream& os) const
{
    os  << "sampledInterface: " << name() ;
}


// ************************************************************************* //
