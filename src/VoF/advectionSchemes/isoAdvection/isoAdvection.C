/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2017 DHI
    Modified code Copyright (C) 2016-2017 OpenCFD Ltd.
    Modified code Copyright (C) 2019 Johan Roenby
    Modified code Copyright (C) 2019-2020 DLR
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

#include "isoAdvection.H"
#include "volFields.H"
#include "interpolationCellPoint.H"
#include "volPointInterpolation.H"
#include "fvcSurfaceIntegrate.H"
#include "fvcGrad.H"
#include "upwind.H"
#include "cellSet.H"
#include "meshTools.H"
#include "OBJstream.H"
#include "syncTools.H"
#include "profiling.H"

#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace advection
{
    defineTypeNameAndDebug(isoAdvection, 0);
    addToRunTimeSelectionTable(advectionSchemes,isoAdvection, components);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::advection::isoAdvection::isoAdvection
(
    volScalarField& alpha1,
    const surfaceScalarField& phi,
    const volVectorField& U
)
:
    advectionSchemes
    (
        typeName,
        alpha1,
        phi,
        U
    ),
    // General data
    mesh_(alpha1.mesh()),
    alpha1In_(alpha1.ref()),
    dVf_
    (
        IOobject
        (
            "dVf_",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimVol, Zero)
    ),
    advectionTime_(0),
    timeIndex_(-1),

    // Tolerances and solution controls
    nAlphaBounds_(modelDict().lookupOrDefault<label>("nAlphaBounds", 3)),
    isoFaceTol_(modelDict().lookupOrDefault<scalar>("isoFaceTol", 1e-8)),
    surfCellTol_(modelDict().lookupOrDefault<scalar>("surfCellTol", 1e-8)),
    writeIsoFacesToFile_(modelDict().lookupOrDefault("writeIsoFaces", false)),

    // Cell cutting data
    surfCells_(label(0.2*mesh_.nCells())),
    advectFace_(alpha1.mesh(), alpha1),
    bsFaces_(label(0.2*mesh_.nBoundaryFaces())),
    bsx0_(bsFaces_.size()),
    bsn0_(bsFaces_.size()),
    bsUn0_(bsFaces_.size()),

    // Parallel run data
    procPatchLabels_(mesh_.boundary().size()),
    surfaceCellFacesOnProcPatches_(0)
{
   cutFaceAdvect::debug = debug;

    // Prepare lists used in parallel runs
    if (Pstream::parRun())
    {
        // Force calculation of required demand driven data (else parallel
        // communication may crash)
        mesh_.cellCentres();
        mesh_.cellVolumes();
        mesh_.faceCentres();
        mesh_.faceAreas();
        mesh_.magSf();
        mesh_.boundaryMesh().patchID();
        mesh_.cellPoints();
        mesh_.cellCells();
        mesh_.cells();

        // Get boundary mesh and resize the list for parallel comms
        setProcessorPatches();
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::advection::isoAdvection::setProcessorPatches()
{
    const polyBoundaryMesh& patches = mesh_.boundaryMesh();
    surfaceCellFacesOnProcPatches_.clear();
    surfaceCellFacesOnProcPatches_.resize(patches.size());

    // Append all processor patch labels to the list
    procPatchLabels_.clear();
    forAll(patches, patchi)
    {
        if
        (
            isA<processorPolyPatch>(patches[patchi])
            && patches[patchi].size() > 0
        )
        {
            procPatchLabels_.append(patchi);
        }
    }
}

void Foam::advection::isoAdvection::extendMarkedCells
(
    bitSet& markedCell
) const
{
    // Mark faces using any marked cell
    bitSet markedFace(mesh_.nFaces());

    for (const label celli : markedCell)
    {
        markedFace.set(mesh_.cells()[celli]);  // set multiple faces
    }

    syncTools::syncFaceList(mesh_, markedFace, orEqOp<unsigned int>());

    // Update cells using any markedFace
    for (label facei = 0; facei < mesh_.nInternalFaces(); ++facei)
    {
        if (markedFace.test(facei))
        {
            markedCell.set(mesh_.faceOwner()[facei]);
            markedCell.set(mesh_.faceNeighbour()[facei]);
        }
    }
    for (label facei = mesh_.nInternalFaces(); facei < mesh_.nFaces(); ++facei)
    {
        if (markedFace.test(facei))
        {
            markedCell.set(mesh_.faceOwner()[facei]);
        }
    }
}

void Foam::advection::isoAdvection::timeIntegratedFlux()
{
    // Get time step
    addProfilingInFunction(geometricVoF);
    const scalar dt = mesh_.time().deltaTValue();

    // Create object for interpolating velocity to isoface centres
    interpolationCellPoint<vector> UInterp(U_);

    // For each downwind face of each surface cell we "isoadvect" to find dVf
    label nSurfaceCells = 0;

    // Clear out the data for re-use and reset list containing information
    // whether cells could possibly need bounding
    clearIsoFaceData();

    // Get necessary references
    const scalarField& phiIn = phi_.primitiveField();
    const scalarField& magSfIn = mesh_.magSf().primitiveField();
    scalarField& dVfIn = dVf_.primitiveFieldRef();

    // Get necessary mesh data
    const cellList& cellFaces = mesh_.cells();
    const labelList& own = mesh_.faceOwner();


    // Storage for isoFace points. Only used if writeIsoFacesToFile_
    DynamicList<List<point>> isoFacePts;
    const DynamicField<label>& interfaceLabels = surf_->interfaceLabels();

    // Loop through cells
    forAll(interfaceLabels, i)
    {
        const label celli = interfaceLabels[i];
        if (mag(surf_->normal()[celli]) != 0)
        {

            // This is a surface cell, increment counter, append and mark cell
            nSurfaceCells++;
            surfCells_.append(celli);

            DebugInfo
                << "\n------------ Cell " << celli << " with alpha1 = "
                << alpha1In_[celli] << " and 1-alpha1 = "
                << 1.0 - alpha1In_[celli] << " ------------"
                << endl;

            // Cell is cut
            const point x0 = surf_->centre()[celli];
            vector n0 = -surf_->normal()[celli];
            n0 /= (mag(n0));

            // Get the speed of the isoface by interpolating velocity and
            // dotting it with isoface unit normal
            const scalar Un0 = UInterp.interpolate(x0, celli) & n0;

            DebugInfo
                << "calcIsoFace gives initial surface: \nx0 = " << x0
                << ", \nn0 = " << n0 << ", \nUn0 = "
                << Un0 << endl;

            // Estimate time integrated flux through each downwind face
            // Note: looping over all cell faces - in reduced-D, some of
            //       these faces will be on empty patches
            const cell& celliFaces = cellFaces[celli];
            forAll(celliFaces, fi)
            {
                const label facei = celliFaces[fi];

                if (mesh_.isInternalFace(facei))
                {
                    bool isDownwindFace = false;

                    if (celli == own[facei])
                    {
                        if (phiIn[facei] >= 0)
                        {
                            isDownwindFace = true;
                        }
                    }
                    else
                    {
                        if (phiIn[facei] < 0)
                        {
                            isDownwindFace = true;
                        }
                    }

                    if (isDownwindFace)
                    {
                        dVfIn[facei] = advectFace_.timeIntegratedFaceFlux
                        (
                            facei,
                            x0,
                            n0,
                            Un0,
                            dt,
                            phiIn[facei],
                            magSfIn[facei]
                        );
                    }

                }
                else
                {
                    bsFaces_.append(facei);
                    bsx0_.append(x0);
                    bsn0_.append(n0);
                    bsUn0_.append(Un0);

                    // Note: we must not check if the face is on the
                    // processor patch here.
                }
            }
        }
    }

    // Get references to boundary fields
    const polyBoundaryMesh& boundaryMesh = mesh_.boundaryMesh();
    const surfaceScalarField::Boundary& phib = phi_.boundaryField();
    const surfaceScalarField::Boundary& magSfb = mesh_.magSf().boundaryField();
    surfaceScalarField::Boundary& dVfb = dVf_.boundaryFieldRef();
    const label nInternalFaces = mesh_.nInternalFaces();

    // Loop through boundary surface faces
    forAll(bsFaces_, i)
    {
        // Get boundary face index (in the global list)
        const label facei = bsFaces_[i];
        const label patchi = boundaryMesh.patchID()[facei - nInternalFaces];
        const label start = boundaryMesh[patchi].start();

        if (phib[patchi].size())
        {
            const label patchFacei = facei - start;
            const scalar phiP = phib[patchi][patchFacei];

            if (phiP >= 0)
            {
                const scalar magSf = magSfb[patchi][patchFacei];

                dVfb[patchi][patchFacei] = advectFace_.timeIntegratedFaceFlux
                (
                    facei,
                    bsx0_[i],
                    bsn0_[i],
                    bsUn0_[i],
                    dt,
                    phiP,
                    magSf
                );

                // Check if the face is on processor patch and append it to
                // the list if necessary
                checkIfOnProcPatch(facei);
            }
        }
    }

    // Synchronize processor patches
    syncProcPatches(dVf_, phi_);

    DebugInfo << "Number of isoAdvector surface cells = "
        << returnReduce(nSurfaceCells, sumOp<label>()) << endl;
}


void Foam::advection::isoAdvection::setDownwindFaces
(
    const label celli,
    DynamicLabelList& downwindFaces
) const
{
    DebugInFunction << endl;

    // Get necessary mesh data and cell information
    const labelList& own = mesh_.faceOwner();
    const cellList& cells = mesh_.cells();
    const cell& c = cells[celli];

    downwindFaces.clear();

    // Check all faces of the cell
    forAll(c, fi)
    {
        // Get face and corresponding flux
        const label facei = c[fi];
        const scalar phi = faceValue(phi_, facei);

        if (own[facei] == celli)
        {
            if (phi >= 0)
            {
                downwindFaces.append(facei);
            }
        }
        else if (phi < 0)
        {
            downwindFaces.append(facei);
        }
    }

    downwindFaces.shrink();
}


Foam::scalar Foam::advection::isoAdvection::netFlux
(
    const surfaceScalarField& dVf,
    const label celli
) const
{
    scalar dV = 0;

    // Get face indices
    const cell& c = mesh_.cells()[celli];

    // Get mesh data
    const labelList& own = mesh_.faceOwner();

    forAll(c, fi)
    {
        const label facei = c[fi];
        const scalar dVff = faceValue(dVf, facei);

        if (own[facei] == celli)
        {
            dV += dVff;
        }
        else
        {
            dV -= dVff;
        }
    }

    return dV;
}


Foam::DynamicList<Foam::label>  Foam::advection::isoAdvection::syncProcPatches
(
    surfaceScalarField& dVf,
    const surfaceScalarField& phi,
    bool returnSyncedFaces
)
{
    DynamicLabelList syncedFaces(0);
    const polyBoundaryMesh& patches = mesh_.boundaryMesh();

    if (Pstream::parRun())
    {
        PstreamBuffers pBufs(Pstream::commsTypes::nonBlocking);

        // Send
        forAll(procPatchLabels_, i)
        {
            const label patchi = procPatchLabels_[i];

            const processorPolyPatch& procPatch =
                refCast<const processorPolyPatch>(patches[patchi]);

            UOPstream toNbr(procPatch.neighbProcNo(), pBufs);
            const scalarField& pFlux = dVf.boundaryField()[patchi];

            const List<label>& surfCellFacesOnProcPatch =
                surfaceCellFacesOnProcPatches_[patchi];

            const UIndirectList<scalar> dVfPatch
            (
                pFlux,
                surfCellFacesOnProcPatch
            );

            toNbr << surfCellFacesOnProcPatch << dVfPatch;
        }

        pBufs.finishedSends();


        // Receive and combine
        forAll(procPatchLabels_, patchLabeli)
        {
            const label patchi = procPatchLabels_[patchLabeli];

            const processorPolyPatch& procPatch =
                refCast<const processorPolyPatch>(patches[patchi]);

            UIPstream fromNeighb(procPatch.neighbProcNo(), pBufs);
            List<label> faceIDs;
            List<scalar> nbrdVfs;

            fromNeighb >> faceIDs >> nbrdVfs;
            if (returnSyncedFaces)
            {
                List<label> syncedFaceI(faceIDs);
                for (label& faceI : syncedFaceI)
                {
                    faceI += procPatch.start();
                }
                syncedFaces.append(syncedFaceI);
            }

            if (debug)
            {
                Pout<< "Received at time = " << mesh_.time().value()
                    << ": surfCellFacesOnProcPatch = " << faceIDs << nl
                    << "Received at time = " << mesh_.time().value()
                    << ": dVfPatch = " << nbrdVfs << endl;
            }

            // Combine fluxes
            scalarField& localFlux = dVf.boundaryFieldRef()[patchi];

            forAll(faceIDs, i)
            {
                const label facei = faceIDs[i];
                localFlux[facei] = - nbrdVfs[i];
                if (debug && mag(localFlux[facei] + nbrdVfs[i]) > ROOTVSMALL)
                {
                    Pout<< "localFlux[facei] = " << localFlux[facei]
                        << " and nbrdVfs[i] = " << nbrdVfs[i]
                        << " for facei = " << facei << endl;
                }
            }
        }

        if (debug)
        {
            // Write out results for checking
            forAll(procPatchLabels_, patchLabeli)
            {
                const label patchi = procPatchLabels_[patchLabeli];
                const scalarField& localFlux = dVf.boundaryField()[patchi];
                Pout<< "time = " << mesh_.time().value() << ": localFlux = "
                    << localFlux << endl;
            }
        }

        // Reinitialising list used for minimal parallel communication
        forAll(surfaceCellFacesOnProcPatches_, patchi)
        {
            surfaceCellFacesOnProcPatches_[patchi].clear();
        }
    }

    return syncedFaces;
}


void Foam::advection::isoAdvection::checkIfOnProcPatch(const label facei)
{
    if (!mesh_.isInternalFace(facei))
    {
        const polyBoundaryMesh& pbm = mesh_.boundaryMesh();
        const label patchi = pbm.patchID()[facei - mesh_.nInternalFaces()];

        if (isA<processorPolyPatch>(pbm[patchi]) && pbm[patchi].size())
        {
            const label patchFacei = pbm[patchi].whichFace(facei);
            surfaceCellFacesOnProcPatches_[patchi].append(patchFacei);
        }
    }
}




void Foam::advection::isoAdvection::applyBruteForceBounding()
{
    addProfilingInFunction(geometricVoF);
    bool alpha1Changed = false;

    scalar snapAlphaTol = modelDict().lookupOrDefault<scalar>("snapTol", 0);
    if (snapAlphaTol > 0)
    {
        alpha1_ =
            alpha1_
           *pos0(alpha1_ - snapAlphaTol)
           *neg0(alpha1_ - (1.0 - snapAlphaTol))
          + pos0(alpha1_ - (1.0 - snapAlphaTol));

        alpha1Changed = true;
    }

    if (modelDict().lookupOrDefault("clip", true))
    {
        alpha1_ = min(scalar(1), max(scalar(0), alpha1_));
        alpha1Changed = true;
    }

    if (alpha1Changed)
    {
        alpha1_.correctBoundaryConditions();
    }
}

// ************************************************************************* //
