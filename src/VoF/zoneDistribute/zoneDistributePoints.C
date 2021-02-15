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

#include "zoneDistributePoints.H"
#include "dummyTransform.H"
#include "emptyPolyPatch.H"
#include "processorPolyPatch.H"
#include "syncTools.H"
#include "wedgePolyPatch.H"

#include "globalPoints.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(zoneDistributePoints, 0);
}


Foam::autoPtr<Foam::indirectPrimitivePatch>
Foam::zoneDistributePoints::nonEmptyWedgePatch() const
{
    const polyBoundaryMesh& patches = mesh_.boundaryMesh();

    label nonWedgeEmpty = 0;

    for (const polyPatch& pp : patches)
    {
        if
        (
            !isA<emptyPolyPatch>(pp) &&
            !isA<wedgePolyPatch>(pp)
        )
        {
            nonWedgeEmpty += pp.size();
        }
    }
    labelList nonWedgeEmptyFaces(nonWedgeEmpty);
    nonWedgeEmpty = 0;

    for (const polyPatch& pp : patches)
    {
        if
        (
            !isA<emptyPolyPatch>(pp) &&
            !isA<wedgePolyPatch>(pp)
        )
        {
            label facei = pp.start();

            forAll(pp, i)
            {
                nonWedgeEmptyFaces[nonWedgeEmpty++] = facei++;
            }
        }
    }

    return autoPtr<indirectPrimitivePatch>::New
    (
        IndirectList<face>
        (
            mesh_.faces(),
            nonWedgeEmptyFaces
        ),
        mesh_.points()
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


Foam::zoneDistributePoints::zoneDistributePoints(const fvMesh& mesh)
:
    MeshObject<fvMesh, Foam::TopologicalMeshObject, zoneDistributePoints>(mesh),
    mesh_(mesh),
    //stencil_(getStencilRef()),
    boundaryPoints_(nonEmptyWedgePatch()().meshPoints()),
    isBoundaryPoint_(mesh.nFaces(),false),
    validBoundaryFace_(mesh.nFaces()-mesh.nInternalFaces(),false)
{

    const polyBoundaryMesh& patches = mesh_.boundaryMesh();

    for (const polyPatch& pp : patches)
    {
        if
        (
            !isA<processorPolyPatch>(pp) &&
            !isA<emptyPolyPatch>(pp) &&
            !isA<wedgePolyPatch>(pp)
        )
        {
            validBoundaryFace_[pp.start() - mesh.nInternalFaces()] = true;
        }
    }
    forAll(boundaryPoints_,i)
    {
        isBoundaryPoint_[boundaryPoints_[i]] = true;
    }
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// void Foam::zoneDistributePoints::updateMesh(const mapPolyMesh& mpm)
// {
//     if(mesh_.topoChanging())
//     {
//         boundaryPoints_ = nonEmptyWedgePatch()().meshPoints();
//         validBoundaryFace_.setSize(mesh_.nFaces()-mesh_.nInternalFaces());
//         validBoundaryFace_ = false;

//         const polyBoundaryMesh& patches = mesh_.boundaryMesh();

//         for (const polyPatch& pp : patches)
//         {
//             if
//             (
//                 !isA<processorPolyPatch>(pp) &&
//                 !isA<emptyPolyPatch>(pp) &&
//                 !isA<wedgePolyPatch>(pp)
//             )
//             {
//                 validBoundaryFace_[pp.start() - mesh_.nInternalFaces()] = true;
//             }
//         }
//         isBoundaryPoint_.setSize(mesh_.nPoints());
//         isBoundaryPoint_ = false;
//         forAll(boundaryPoints_,i)
//         {
//             isBoundaryPoint_[boundaryPoints_[i]] = true;
//         }
//     }

// }

// ************************************************************************* //
