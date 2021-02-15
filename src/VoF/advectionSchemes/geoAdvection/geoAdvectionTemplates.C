/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2017 DHI
    Modified code Copyright (C) 2016-2017 OpenCFD Ltd.
    Modified code Copyright (C) 2019 DLR
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

#include "geoAdvection.H"
#include "fvcSurfaceIntegrate.H"
#include "upwind.H"
#include "interpolationCellPoint.H"

// ************************************************************************* //

template<typename Type>
Type Foam::advection::geoAdvection::faceValue
(
    const GeometricField<Type, fvsPatchField, surfaceMesh>& f,
    const label facei
) const
{
    if (mesh_.isInternalFace(facei))
    {
        return f.primitiveField()[facei];
    }
    else
    {
        const polyBoundaryMesh& pbm = mesh_.boundaryMesh();

        // Boundary face. Find out which face of which patch
        const label patchi = pbm.patchID()[facei - mesh_.nInternalFaces()];

        if (patchi < 0 || patchi >= pbm.size())
        {
            FatalErrorInFunction
                << "Cannot find patch for face " << facei
                << abort(FatalError);
        }

        // Handle empty patches
        const polyPatch& pp = pbm[patchi];
        if (isA<emptyPolyPatch>(pp) || pp.empty())
        {
            return pTraits<Type>::zero;
        }

        const label patchFacei = pp.whichFace(facei);
        return f.boundaryField()[patchi][patchFacei];
    }
}


template<typename Type>
void Foam::advection::geoAdvection::setFaceValue
(
    GeometricField<Type, fvsPatchField, surfaceMesh>& f,
    const label facei,
    const Type& value
) const
{
    if (mesh_.isInternalFace(facei))
    {
        f.primitiveFieldRef()[facei] = value;
    }
    else
    {
        const polyBoundaryMesh& pbm = mesh_.boundaryMesh();

        // Boundary face. Find out which face of which patch
        const label patchi = pbm.patchID()[facei - mesh_.nInternalFaces()];

        if (patchi < 0 || patchi >= pbm.size())
        {
            FatalErrorInFunction
                << "Cannot find patch for face " << facei
                << abort(FatalError);
        }

        // Handle empty patches
        const polyPatch& pp = pbm[patchi];
        if (isA<emptyPolyPatch>(pp) || pp.empty())
        {
            return;
        }

        const label patchFacei = pp.whichFace(facei);

        f.boundaryFieldRef()[patchi][patchFacei] = value;
    }
}

template < class RdeltaTType,class SpType, class SuType >
void Foam::advection::geoAdvection::limitFluxes
(
    const RdeltaTType& rDeltaT,
    const SpType& Sp,
    const SuType& Su
)
{
    DebugInFunction << endl;

    const scalar aTol = 1.0e-12;          // Note: tolerances
    scalar maxAlphaMinus1 = gMax(alpha1_) - 1;      // max(alphaNew - 1);
    scalar minAlpha = gMin(alpha1_);           // min(alphaNew);
    const label nOvershoots = 20;         // sum(pos0(alphaNew - 1 - aTol));

    const labelList& owner = mesh_.faceOwner();
    const labelList& neighbour = mesh_.faceNeighbour();

    DebugInfo << "geoAdvection: Before conservative bounding: min(alpha) = "
        << minAlpha << ", max(alpha) = 1 + " << maxAlphaMinus1 << endl;

    surfaceScalarField afcorrectionValues("afcorrectionValues", alphaPhi_*0.0);


    // Loop number of bounding steps
    for (label n = 0; n < nAlphaBounds_; n++)
    {
        if (maxAlphaMinus1 > aTol || minAlpha < -aTol) // Note: tolerances
        {
            DebugInfo << "boundAlpha... " << endl;

            DynamicList<label> correctedFaces(3*nOvershoots);
            afcorrectionValues = dimensionedScalar("0",dimVolume/dimTime,0.0);
            boundFlux(alpha1In_, afcorrectionValues, correctedFaces,rDeltaT,Sp,Su);

            correctedFaces.append
            (
                syncProcPatches(afcorrectionValues, phi_,true)
            );

            labelHashSet alreadyUpdated;
            forAll(correctedFaces, fi)
            {
                label facei = correctedFaces[fi];
                if(alreadyUpdated.insert(facei))
                {
                    checkIfOnProcPatch(facei);
                    const label own = owner[facei];
                    const scalar dt = deltaT(rDeltaT,own);

                    alpha1_[own] +=
                        -faceValue(afcorrectionValues, facei)*dt/mesh_.V()[own];
                    if(mesh_.isInternalFace(facei))
                    {
                        const label nei = neighbour[facei];
                        alpha1_[nei] -=
                            - faceValue(afcorrectionValues, facei)*dt
                            / mesh_.V()[nei];
                    }

                    // Change to treat boundaries consistently
                    scalar corralphaf =
                        faceValue(alphaPhi_, facei)
                      + faceValue(afcorrectionValues, facei);

                    setFaceValue(alphaPhi_, facei, corralphaf);
                }
            }
            syncProcPatches(alphaPhi_, phi_);
        }
        else
        {
            break;
        }

        maxAlphaMinus1 = gMax(alpha1_) - 1;     // max(alphaNew - 1);
        minAlpha = gMin(alpha1_);               // min(alphaNew);

        if (debug)
        {
            // Check if still unbounded
            //scalarField alphaNew(alpha1In_ - fvc::surfaceIntegrate(alphaPhi_)());
            label maxAlphaMinus1 = max(alpha1_.primitiveField() - 1);
            scalar minAlpha = min(alpha1_.primitiveField());
            label nUndershoots = sum(neg0(alpha1_.primitiveField() + aTol));
            label nOvershoots = sum(pos0(alpha1_.primitiveField() - 1 - aTol));

            Info<< "After bounding number " << n + 1 << " of time "
                << mesh_.time().value() << ":" << endl;
            Info<< "nOvershoots = " << nOvershoots << " with max(alpha1_-1) = "
                << maxAlphaMinus1 << " and nUndershoots = " << nUndershoots
                << " with min(alpha1_) = " << minAlpha << endl;
        }
    }

    alpha1_.correctBoundaryConditions();

}


template < class RdeltaTType,class SpType, class SuType >
void Foam::advection::geoAdvection::boundFlux
(
    const scalarField& alpha1,
    surfaceScalarField& afcorrectionValues,
    DynamicList<label>& correctedFaces,
    const RdeltaTType& rDeltaT,
    const SpType& Sp,
    const SuType& Su
)
{
    DebugInFunction << endl;

    correctedFaces.clear();
    scalar aTol = 10*SMALL; // Note: tolerances

    const scalarField& meshV = mesh_.cellVolumes();

    DynamicList<label> downwindFaces(10);
    DynamicList<label> facesToPassFluidThrough(downwindFaces.size());
    DynamicList<scalar> dVfmax(downwindFaces.size());
    DynamicList<scalar> phi(downwindFaces.size());
    const volScalarField& alphaOld = alpha1_.oldTime();

    // Loop through alpha cell centred field
    forAll(alpha1, celli)
    {
        if (alpha1[celli] < -aTol || alpha1[celli] > 1 + aTol)
        {
            const scalar dt = deltaT(rDeltaT,celli);
            const scalar rDt = 1/dt;
            const scalar Vi = meshV[celli];
            scalar alphaOvershoot = pos0(alpha1[celli]-1)*(alpha1[celli]-1)
                + neg0(alpha1[celli])*alpha1[celli];
            scalar fluidToPassOn = alphaOvershoot*Vi;
            label nFacesToPassFluidThrough = 1;

            bool firstLoop = true;

            // First try to pass surplus fluid on to neighbour cells that are
            // not filled and to which dVf < phi*dt
            while (mag(alphaOvershoot) > aTol && nFacesToPassFluidThrough > 0)
            {
                DebugInfo
                    << "\n\nBounding cell " << celli
                    << " with alpha overshooting " << alphaOvershoot
                    << endl;

                facesToPassFluidThrough.clear();
                dVfmax.clear();
                phi.clear();

                // Find potential neighbour cells to pass surplus phase to
                setDownwindFaces(celli, downwindFaces);

                scalar dVftot = 0;
                nFacesToPassFluidThrough = 0;

                forAll(downwindFaces, fi)
                {
                    const label facei = downwindFaces[fi];
                    const scalar phif = faceValue(phi_, facei);

                    const scalar dVff =
                        faceValue(alphaPhi_, facei)*dt
                      + faceValue(afcorrectionValues, facei)*dt;

                    const scalar maxExtraFaceFluidTrans =
                        mag(pos0(fluidToPassOn)*phif*dt - dVff);

                    // dVf has same sign as phi and so if phi>0 we have
                    // mag(phi_[facei]*dt) - mag(dVf[facei]) = phi_[facei]*dt
                    // - dVf[facei]
                    // If phi < 0 we have mag(phi_[facei]*dt) -
                    // mag(dVf[facei]) = -phi_[facei]*dt - (-dVf[facei]) > 0
                    // since mag(dVf) < phi*dt
                    DebugInfo
                        << "downwindFace " << facei
                        << " has maxExtraFaceFluidTrans = "
                        << maxExtraFaceFluidTrans << endl;

                    if (maxExtraFaceFluidTrans/Vi > aTol)
                    {
//                    if (maxExtraFaceFluidTrans/Vi > aTol &&
//                    mag(dVfIn[facei])/Vi > aTol) //Last condition may be
//                    important because without this we will flux through uncut
//                    downwind faces
                        facesToPassFluidThrough.append(facei);
                        phi.append(phif);
                        dVfmax.append(maxExtraFaceFluidTrans);
                        dVftot += mag(phif*dt);
                    }
                }

                DebugInfo
                    << "\nfacesToPassFluidThrough: "
                    << facesToPassFluidThrough << ", dVftot = "
                    << dVftot << " m3 corresponding to dalpha = "
                    << dVftot/Vi << endl;

                forAll(facesToPassFluidThrough, fi)
                {
                    const label facei = facesToPassFluidThrough[fi];
                    scalar fluidToPassThroughFace =
                        mag(fluidToPassOn)*mag(phi[fi]*dt)/dVftot;

                    nFacesToPassFluidThrough +=
                        pos0(dVfmax[fi] - fluidToPassThroughFace);

                    fluidToPassThroughFace =
                        min(fluidToPassThroughFace, dVfmax[fi]);

                    scalar dVff = faceValue(afcorrectionValues, facei)*dt;

                    dVff +=
                        sign(phi[fi])*fluidToPassThroughFace*sign(fluidToPassOn);

                    setFaceValue(afcorrectionValues, facei, dVff/dt);

                    if (firstLoop)
                    {
                        checkIfOnProcPatch(facei);
                        correctedFaces.append(facei);
                    }
                }

                firstLoop = false;

                scalar alpha1New =
                (
                    alphaOld[celli]*rDt  + Su[celli]
                  - netFlux(alphaPhi_, celli)/Vi
                  - netFlux(afcorrectionValues, celli)/Vi
                )
                /
                (rDt - Sp[celli]);

                alphaOvershoot =
                    pos0(alpha1New-1)*(alpha1New-1)
                  + neg0(alpha1New)*alpha1New;

                fluidToPassOn = alphaOvershoot*Vi;

                DebugInfo
                    << "\nNew alpha for cell " << celli << ": "
                    << alpha1New << endl;
            }
        }
    }

    DebugInfo << "correctedFaces = " << correctedFaces << endl;
}

template < class SpType, class SuType >
void Foam::advection::geoAdvection::advect
(
    const SpType& Sp,
    const SuType& Su
)
{
    if (fv::localEulerDdt::enabled(mesh_))
    {
        const volScalarField& rDeltaT = fv::localEulerDdt::localRDeltaT(mesh_);
        advect(rDeltaT,Sp,Su);
    }
    else
    {
        scalar rDeltaT = 1/mesh_.time().deltaTValue();
        advect(rDeltaT,Sp,Su);
    }

}


template < class RdeltaTType >
void Foam::advection::geoAdvection::timeIntegratedFlux
(
    const RdeltaTType& rDeltaT
)
{
    // Get time step

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
    scalarField& alphaPhiIn = alphaPhi_.primitiveFieldRef();

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
            const scalar dt = deltaT(rDeltaT,celli);
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
                        alphaPhiIn[facei] = advectFace_.timeIntegratedFaceFlux
                        (
                            facei,
                            x0,
                            n0,
                            Un0,
                            dt,
                            phiIn[facei],
                            magSfIn[facei]
                        )/dt;
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
    surfaceScalarField::Boundary& alphaPhib = alphaPhi_.boundaryFieldRef();
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
            const label celli = mesh_.faceOwner()[facei];
            const scalar dt = deltaT(rDeltaT,celli);
            const scalar phiP = phib[patchi][patchFacei];

            if (phiP >= 0)
            {
                const scalar magSf = magSfb[patchi][patchFacei];

                alphaPhib[patchi][patchFacei] = advectFace_.timeIntegratedFaceFlux
                (
                    facei,
                    bsx0_[i],
                    bsn0_[i],
                    bsUn0_[i],
                    dt,
                    phiP,
                    magSf
                )/dt;

                // Check if the face is on processor patch and append it to
                // the list if necessary
                checkIfOnProcPatch(facei);
            }
        }
    }

    // Synchronize processor patches
    syncProcPatches(alphaPhi_, phi_);

    DebugInfo << "Number of isoAdvector surface cells = "
        << returnReduce(nSurfaceCells, sumOp<label>()) << endl;
}


template < class RdeltaTType,class SpType, class SuType >
void Foam::advection::geoAdvection::advect
(
    const RdeltaTType& rDeltaT,
    const SpType& Sp,
    const SuType& Su
)
{
    DebugInFunction << endl;

    if (mesh_.topoChanging())
    {
        setProcessorPatches();
    }

    scalar advectionStartTime = mesh_.time().elapsedCpuTime();

    // reconstruct the interface
    surf_->reconstruct();

    // Initialising alphaPhi with upwind values
    // i.e. phi[facei]*alpha1[upwindCell[facei]]
    alphaPhi_ = upwind<scalar>(mesh_, phi_).flux(alpha1_);


    // Do the geoAdvection on surface cells
    timeIntegratedFlux(rDeltaT);

    // Adjust alpha for mesh motion
    if (mesh_.moving())
    {
        alpha1In_ *= (mesh_.Vsc0()/mesh_.Vsc());
    }

    // Advect the free surface
    alpha1_.primitiveFieldRef() =
    (
        alpha1_.oldTime().primitiveField()*rDeltaT
      + Su.field()
      - fvc::surfaceIntegrate(alphaPhi_)().primitiveField()
    )/(rDeltaT - Sp.field());

    alpha1_.correctBoundaryConditions();

    // Adjust alphaPhi for unbounded cells
    limitFluxes
    (
        rDeltaT,
        Sp,
        Su
    );

    scalar maxAlphaMinus1 = gMax(alpha1In_) - 1;
    scalar minAlpha = gMin(alpha1In_);

    Info<< "geoAdvection: After conservative bounding: min(alpha) = "
        << minAlpha << ", max(alpha) = 1 + " << maxAlphaMinus1 << endl;

    // Apply non-conservative bounding mechanisms (clipping and snapping)
    // Note: We should be able to write out alpha before this is done!
    applyBruteForceBounding();

    advectionTime_ += (mesh_.time().elapsedCpuTime() - advectionStartTime);
    DebugInfo
        << "geoAdvection: time consumption = "
        << label(100*advectionTime_/(mesh_.time().elapsedCpuTime() + SMALL))
        << "%" << endl;

}
// ************************************************************************* //
