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

//#include "dummyTransform.H"
#include "emptyPolyPatch.H"
#include "markInterfaceRegion.H"
#include "processorPolyPatch.H"
#include "syncTools.H"
#include "wedgePolyPatch.H"
#include "indexedOctree.H"
#include "treeDataPoint.H"
#include "alphaContactAngleTwoPhaseFvPatchScalarField.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::autoPtr<Foam::indirectPrimitivePatch>
Foam::markInterfaceRegion::coupledFacesPatch() const
{
    const polyBoundaryMesh& patches = mesh_.boundaryMesh();

    label nCoupled = 0;

    for (const polyPatch& pp : patches)
    {
        if (isA<coupledPolyPatch>(pp))
        {
            nCoupled += pp.size();
        }
    }
    labelList nCoupledFaces(nCoupled);
    nCoupled = 0;

    for (const polyPatch& pp : patches)
    {
        if (isA<coupledPolyPatch>(pp))
        {
            label facei = pp.start();

            forAll(pp, i)
            {
                nCoupledFaces[nCoupled++] = facei++;
            }
        }
    }

    return autoPtr<indirectPrimitivePatch>::New
    (
        IndirectList<face>
        (
            mesh_.faces(),
            nCoupledFaces
        ),
        mesh_.points()
    );
}

void Foam::markInterfaceRegion::markCellsNearSurf
(
    const boolList& interfaceCells,
    const label neiRingLevel,
    boolList& nextToInterface,
    volScalarField& cellDistLevel
)
{
    // performance might be improved by increasing the saving last iteations
    // cells in a Map and loop over the map
    if (mesh_.topoChanging())
    {
        // Introduced resizing to cope with changing meshes
        if (nextToInterface.size() != mesh_.nCells())
        {
            nextToInterface.resize(mesh_.nCells());
        }
        coupledBoundaryPoints_ = coupledFacesPatch()().meshPoints();
    }

    const labelListList& pCells = mesh_.cellPoints();
    const labelListList& cPoints = mesh_.pointCells();

    boolList alreadyMarkedPoint(mesh_.nPoints(),false);
    nextToInterface = false;

    // do coupled face first
    Map<bool> syncMap;

    for(int level = 0;level<=neiRingLevel;level++)
    {
        // parallel
        if (level > 0)
        {
            forAll(coupledBoundaryPoints_,i)
            {
                const label& pi = coupledBoundaryPoints_[i];
                forAll(mesh_.pointCells()[pi], j)
                {
                    const label celli = cPoints[pi][j];
                    if (cellDistLevel[celli] == level-1)
                    {
                        syncMap.insert(pi, true);
                        break;
                    }
                }
            }

            syncTools::syncPointMap(mesh_, syncMap, orEqOp<bool>());

            // mark parallel points first
            forAllConstIter(Map<bool>, syncMap, iter)
            {
                const label pi = iter.key();

                if (!alreadyMarkedPoint[pi])
                {
                    // loop over all cells attached to the point
                    forAll(cPoints[pi],j)
                    {
                        const label& pCelli = cPoints[pi][j];
                        if (cellDistLevel[pCelli] == -1)
                        {
                            cellDistLevel[pCelli] = level;
                            nextToInterface[pCelli] = true;
                        }

                    }
                }

                alreadyMarkedPoint[pi] = true;
            }
        }


        forAll(cellDistLevel, celli)
        {
            if (level == 0)
            {
                if (interfaceCells[celli])
                {
                    cellDistLevel[celli] = 0;
                    nextToInterface[celli] = true;
                }
                else
                {
                    cellDistLevel[celli] = -1;
                }
            }
            else
            {
                if (cellDistLevel[celli] == level-1)
                {
                    forAll(pCells[celli], i)
                    {
                        const label& pI = pCells[celli][i];
                        // interfacePoint_[pI] = true;
                        if (!alreadyMarkedPoint[pI])
                        {
                            // loop over all cells attached to the point
                            forAll(cPoints[pI],j)
                            {
                                const label& pCelli = cPoints[pI][j];
                                if (cellDistLevel[pCelli] == -1)
                                {
                                    cellDistLevel[pCelli] = level;
                                    nextToInterface[pCelli] = true;
                                }

                            }
                        }
                        alreadyMarkedPoint[pI] = true;

                    }
                }
            }
        }

    }

}


void Foam::markInterfaceRegion::markCellsNearSurf
(
    const boolList& interfaceCells,
    const label neiRingLevel,
    boolList& nextToInterface
)
{
    // performance might be improved by increasing the saving last iteations
    // cells in a Map and loop over the map
    if (mesh_.topoChanging())
    {
        // Introduced resizing to cope with changing meshes
        if (nextToInterface.size() != mesh_.nCells())
        {
            nextToInterface.resize(mesh_.nCells());
        }
        coupledBoundaryPoints_ = coupledFacesPatch()().meshPoints();
    }
    labelField cellDistLevel(mesh_.nCells());

    const labelListList& pCells = mesh_.cellPoints();
    const labelListList& cPoints = mesh_.pointCells();

    boolList alreadyMarkedPoint(mesh_.nPoints(),false);
    nextToInterface = false;

    // do coupled face first
    Map<bool> syncMap;

    for(int level = 0;level<=neiRingLevel;level++)
    {
        // parallel
        if (level > 0)
        {
            forAll(coupledBoundaryPoints_,i)
            {
                const label& pi = coupledBoundaryPoints_[i];
                forAll(mesh_.pointCells()[pi], j)
                {
                    const label celli = cPoints[pi][j];
                    if (cellDistLevel[celli] == level-1)
                    {
                        syncMap.insert(pi, true);
                        break;
                    }
                }
            }

            syncTools::syncPointMap(mesh_, syncMap, orEqOp<bool>());

            // mark parallel points first
            forAllConstIter(Map<bool>, syncMap, iter)
            {
                const label pi = iter.key();

                if (!alreadyMarkedPoint[pi])
                {
                    // loop over all cells attached to the point
                    forAll(cPoints[pi],j)
                    {
                        const label& pCelli = cPoints[pi][j];
                        if (cellDistLevel[pCelli] == -1)
                        {
                            cellDistLevel[pCelli] = level;
                            nextToInterface[pCelli] = true;
                        }

                    }
                }

                alreadyMarkedPoint[pi] = true;
            }
        }


        forAll(cellDistLevel, celli)
        {
            if (level == 0)
            {
                if (interfaceCells[celli])
                {
                    cellDistLevel[celli] = 0;
                    nextToInterface[celli] = true;
                }
                else
                {
                    cellDistLevel[celli] = -1;
                }
            }
            else
            {
                if (cellDistLevel[celli] == level-1)
                {
                    forAll(pCells[celli], i)
                    {
                        const label& pI = pCells[celli][i];
                        // interfacePoint_[pI] = true;
                        if (!alreadyMarkedPoint[pI])
                        {
                             // loop over all cells attached to the point
                            forAll(cPoints[pI],j)
                            {
                                const label& pCelli = cPoints[pI][j];
                                if (cellDistLevel[pCelli] == -1)
                                {
                                    cellDistLevel[pCelli] = level;
                                    nextToInterface[pCelli] = true;
                                }

                            }
                        }
                        alreadyMarkedPoint[pI] = true;

                    }
                }
            }
        }

    }

}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::markInterfaceRegion::markInterfaceRegion
(
    const fvMesh& mesh
)
:
    mesh_(mesh),
    coupledBoundaryPoints_(coupledFacesPatch()().meshPoints())
{
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


// ************************************************************************* //
