/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2019-2019 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
                            | Copyright (C) 2019-2019 DLR
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
#include "dummyTransform.H"
#include "emptyPolyPatch.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(zoneCPCStencil, 0);
}


Foam::Map<bool> Foam::zoneCPCStencil::syncCouledBoundaryPoints
(
    const boolList& zone,
    const labelList& boundaryPoints
) const
{
    const labelListList& pCells = mesh_.pointCells();

    Map<bool> syncPoints;

    forAll(boundaryPoints, i)
    {
        label pointi = boundaryPoints[i];

        bool updatePoint = false;

        // check if point need to be updated
        forAll(pCells[pointi],j)
        {
            const label& celli = pCells[pointi][j];
            if(zone[celli])
            {
                updatePoint = true;
                break;
            }

        }

        if(updatePoint)
        {
            syncPoints.insert(pointi,true);
        }
    }

    // sync syncPoints
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

    forAll(boundaryPoints, i)
    {
        label pointi = boundaryPoints[i];

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
        Foam::dummyTransform()      // dummy transformation
    );
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::zoneCPCStencil::zoneCPCStencil(const fvMesh& mesh)
:
    MeshObject<fvMesh, Foam::TopologicalMeshObject, zoneCPCStencil>(mesh),
    zoneCellStencils(mesh),
    nonEmptyBoundaryPoints_(nonEmptyFacesPatch()().meshPoints()),
    uptodate_(mesh.nCells(),false)
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
    Map<bool> syncPoints  = syncCouledBoundaryPoints
    (
        zone,
        nonEmptyBoundaryPoints_
    );

    labelList boundaryPoints = syncPoints.toc();

    Map<labelList> neiGlobal;
    calcPointBoundaryData
    (
        zone,
        isValidBFace_,
        boundaryPoints,
        neiGlobal
    );

    // add boundary Points first

    forAll(boundaryPoints, i)
    {
        label pointi = boundaryPoints[i];

        const labelList& pGlobals = neiGlobal[pointi];

        // Distribute to all pointCells
        const labelList& pCells = mesh_.pointCells(pointi);

        forAll(pCells, j)
        {
            label celli = pCells[j];

            // Insert pGlobals into globalCellCells
            if(zone[celli] && !uptodate_[celli])
            {
                merge
                (
                    globalNumbering().toGlobal(celli),
                    pGlobals,
                    globalCellCells[celli]
                );

                forAll(globalCellCells[celli],idx)
                {
                    const label& gblIdx = globalCellCells[celli][idx];
                    if(!globalNumbering().isLocal(gblIdx))
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

        if(zone[celli] && !uptodate_[celli])
        {
            forAll(cPoints[celli],j)
            {
                const label& pointi = cPoints[celli][j];

                labelList pCells = mesh_.pointCells(pointi);

                forAll(pCells, cpi)
                {
                    label neiCelli = pCells[cpi];
                    pCells[cpi] = globalNumbering().toGlobal(neiCelli);
                }

                if(!uptodate_[celli])
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

// void Foam::zoneCPCStencil::updateMesh(const mapPolyMesh& mpm)
// {
//     if(mesh_.topoChanging())
//     { 
//         // resize map and globalIndex
//         zoneCellStencils::updateMesh(mpm);
        
//         nonEmptyBoundaryPoints_ = nonEmptyFacesPatch()().meshPoints();
//         uptodate_.resize(mesh_.nCells());
//         uptodate_ = false;
//         validBoundaryFaces(isValidBFace_);
//     } 

// }

// ************************************************************************* //
