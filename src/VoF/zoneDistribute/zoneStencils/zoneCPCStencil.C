/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020 DLR
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

#include "zoneCPCStencil.H"
#include "syncTools.H"
#include "ListOps.H"
#include "dummyTransform.H"
#include "emptyPolyPatch.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(zoneCPCStencil, 0);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::Map<bool> Foam::zoneCPCStencil::syncCoupledBoundaryPoints
(
    const boolList& zone,
    const labelList& boundaryPoints
) const
{
    const labelListList& pCells = mesh_.pointCells();

    Map<bool> syncPoints;

    for (const label pointi : boundaryPoints)
    {
        bool updatePoint = false;

        // Check if point needs to be updated
        for (const label celli : pCells[pointi])
        {
            if (zone[celli])
            {
                updatePoint = true;
                break;
            }
        }

        if (updatePoint)
        {
            syncPoints.insert(pointi, true);
        }
    }

    syncTools::syncPointMap
    (
        mesh_,
        syncPoints,
        orEqOp<bool>(),
        Foam::dummyTransform()
    );

    return syncPoints;
}


void Foam::zoneCPCStencil::calcPointBoundaryData
(
    const boolList& zone,
    const boolList& isValidBFace,
    const labelList& boundaryPoints,
    Map<labelList>& neiGlobal
) const
{
    neiGlobal.resize(2*boundaryPoints.size());

    labelHashSet pointGlobals;

    for (const label pointi : boundaryPoints)
    {
        neiGlobal.insert
        (
            pointi,
            calcFaceCells
            (
                isValidBFace,
                mesh_.pointFaces()[pointi],
                pointGlobals
            )
        );
    }

    syncTools::syncPointMap
    (
        mesh_,
        neiGlobal,
        unionEqOp(),
        Foam::dummyTransform()
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::zoneCPCStencil::zoneCPCStencil(const fvMesh& mesh)
:
    MeshObject<fvMesh, Foam::TopologicalMeshObject, zoneCPCStencil>(mesh),
    zoneCellStencils(mesh),
    nonEmptyBoundaryPoints_(nonEmptyFacesPatch()().meshPoints()),
    uptodate_(mesh.nCells(), false)
{
    // Mark boundary faces to be included in stencil (i.e. not coupled or empty)
    validBoundaryFaces(isValidBFace_);
}

Foam::zoneCPCStencil& Foam::zoneCPCStencil::New(const fvMesh& mesh)
{
    bool found = mesh.thisDb().foundObject<zoneCPCStencil>
    (
        zoneCPCStencil::typeName
    );
    zoneCPCStencil* ptr = nullptr;

    if(found)
    {
        ptr = mesh.thisDb().getObjectPtr<zoneCPCStencil>
        (
            zoneCPCStencil::typeName
        );

        return *ptr;
    }

    zoneCPCStencil* objectPtr = new zoneCPCStencil(mesh);

    regIOobject::store(static_cast<zoneCPCStencil*>(objectPtr));

    return *objectPtr;
}


void Foam::zoneCPCStencil::calculateStencil
(
    const boolList& zone,
    labelListList& globalCellCells
)
{
    // Swap pointCells for coupled points
    Map<bool> syncPoints = syncCoupledBoundaryPoints
    (
        zone,
        nonEmptyBoundaryPoints_
    );

    labelList boundaryPoints(syncPoints.toc());

    Map<labelList> neiGlobal;
    calcPointBoundaryData
    (
        zone,
        isValidBFace_,
        boundaryPoints,
        neiGlobal
    );

    // add boundary Points first

    for (const label pointi : boundaryPoints)
    {
        const labelList& pGlobals = neiGlobal[pointi];

        // Distribute to all pointCells
        const labelList& pCells = mesh_.pointCells(pointi);

        for (const label celli : pCells)
        {
            // Insert pGlobals into globalCellCells
            if (zone[celli] && !uptodate_[celli])
            {
                merge
                (
                    globalNumbering().toGlobal(celli),
                    pGlobals,
                    globalCellCells[celli]
                );

                for (const label gblIdx : globalCellCells[celli])
                {
                    if (!globalNumbering().isLocal(gblIdx))
                    {
                        needComm_.insert(celli);
                    }
                }
            }
        }
    }


    neiGlobal.clear();

    // Do remaining points cells
    const labelListList& cPoints = mesh_.cellPoints();

    forAll(zone,celli)
    {
        if (zone[celli] && !uptodate_[celli])
        {
            for (const label pointi : cPoints[celli])
            {
                labelList pCells = mesh_.pointCells(pointi);

                for (label& neiCelli : pCells)
                {
                    neiCelli = globalNumbering().toGlobal(neiCelli);
                }

                if (!uptodate_[celli])
                {
                    merge
                    (
                        globalNumbering().toGlobal(celli),
                        pCells,
                        globalCellCells[celli]
                    );
                }
            }

            uptodate_[celli] = true;
        }
    }
}


// ************************************************************************* //
