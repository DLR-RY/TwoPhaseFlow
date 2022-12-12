/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2017 DHI
    Modified code Copyright (C) 2016-2017 OpenCFD Ltd.
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
#include "fvcSurfaceIntegrate.H"
#include "upwind.H"

// ************************************************************************* //

template<typename Type>
Type Foam::advection::isoAdvection::faceValue
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
void Foam::advection::isoAdvection::setFaceValue
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


template<class SpType, class SuType>
void Foam::advection::isoAdvection::limitFluxes
(
    const SpType& Sp,
    const SuType& Su
)
{
    addProfilingInFunction(geometricVoF);
    DebugInFunction << endl;

    const scalar aTol = 1.0e-12;          // Note: tolerances
    scalar maxAlphaMinus1 = gMax(alpha1_) - 1;      // max(alphaNew - 1);
    scalar minAlpha = gMin(alpha1_);           // min(alphaNew);
    const label nOvershoots = 20;         // sum(pos0(alphaNew - 1 - aTol));

    const labelList& owner = mesh_.faceOwner();
    const labelList& neighbour = mesh_.faceNeighbour();

    DebugInfo
        << "isoAdvection: Before conservative bounding: min(alpha) = "
        << minAlpha << ", max(alpha) = 1 + " << maxAlphaMinus1 << endl;

    surfaceScalarField dVfcorrectionValues("dVfcorrectionValues", dVf_*0.0);

    bitSet needBounding(mesh_.nCells(),false);
    needBounding.set(surfCells_);

    extendMarkedCells(needBounding);

    // Loop number of bounding steps
    for (label n = 0; n < nAlphaBounds_; n++)
    {
        if (maxAlphaMinus1 > aTol || minAlpha < -aTol) // Note: tolerances
        {
            DebugInfo << "boundAlpha... " << endl;

            DynamicList<label> correctedFaces(3*nOvershoots);
            dVfcorrectionValues = dimensionedScalar("0",dimVolume,0.0);
            boundFlux(needBounding, dVfcorrectionValues, correctedFaces,Sp,Su);

            correctedFaces.append
            (
                syncProcPatches(dVfcorrectionValues, phi_,true)
            );

            labelHashSet alreadyUpdated;
            forAll(correctedFaces, fi)
            {
                label facei = correctedFaces[fi];
                if (alreadyUpdated.insert(facei))
                {
                    checkIfOnProcPatch(facei);
                    const label own = owner[facei];

                    alpha1_[own] +=
                        -faceValue(dVfcorrectionValues, facei)/mesh_.V()[own];
                    if (mesh_.isInternalFace(facei))
                    {
                        const label nei = neighbour[facei];
                        alpha1_[nei] -=
                            -faceValue(dVfcorrectionValues, facei)/mesh_.V()[nei];
                    }

                    // Change to treat boundaries consistently
                    scalar corrVf =
                        faceValue(dVf_, facei)
                      + faceValue(dVfcorrectionValues, facei);

                    setFaceValue(dVf_, facei, corrVf);
                }
            }
            syncProcPatches(dVf_, phi_);
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
            //scalarField alphaNew(alpha1In_ - fvc::surfaceIntegrate(dVf_)());
            label maxAlphaMinus1 = max(alpha1_.primitiveField() - 1);
            scalar minAlpha = min(alpha1_.primitiveField());
            label nUndershoots = sum(neg0(alpha1_.primitiveField() + aTol));
            label nOvershoots = sum(pos0(alpha1_.primitiveField() - 1 - aTol));

            Info<< "After bounding number " << n + 1 << " of time "
                << mesh_.time().value() << ":" << nl
                << "nOvershoots = " << nOvershoots << " with max(alpha1_-1) = "
                << maxAlphaMinus1 << " and nUndershoots = " << nUndershoots
                << " with min(alpha1_) = " << minAlpha << endl;
        }
    }

    alpha1_.correctBoundaryConditions();

}


template<class SpType, class SuType>
void Foam::advection::isoAdvection::boundFlux
(
    const bitSet& nextToInterface,
    surfaceScalarField& dVfcorrectionValues,
    DynamicList<label>& correctedFaces,
    const SpType& Sp,
    const SuType& Su
)
{
    addProfilingInFunction(geometricVoF);
    DebugInFunction << endl;
    scalar rDeltaT = 1/mesh_.time().deltaTValue();

    correctedFaces.clear();
    scalar aTol = 10*SMALL; // Note: tolerances

    const scalarField& meshV = mesh_.cellVolumes();
    const scalar dt = mesh_.time().deltaTValue();

    DynamicList<label> downwindFaces(10);
    DynamicList<label> facesToPassFluidThrough(downwindFaces.size());
    DynamicList<scalar> dVfmax(downwindFaces.size());
    DynamicList<scalar> phi(downwindFaces.size());
    const volScalarField& alphaOld = alpha1_.oldTime();

    // Loop through alpha cell centred field
    for(label celli: nextToInterface)
    {
        if (alpha1_[celli] < -aTol || alpha1_[celli] > 1 + aTol)
        {
            const scalar Vi = meshV[celli];
            scalar alphaOvershoot = pos0(alpha1_[celli]-1)*(alpha1_[celli]-1)
                + neg0(alpha1_[celli])*alpha1_[celli];
            scalar fluidToPassOn = alphaOvershoot*Vi;
            label nFacesToPassFluidThrough = 1;

            bool firstLoop = true;

            // First try to pass surplus fluid on to neighbour cells that are
            // not filled and to which dVf < phi*dt
            for (label iter=0; iter<10; iter++)
            {
                if(mag(alphaOvershoot) < aTol || nFacesToPassFluidThrough == 0)
                {
                    break;
                }

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
                        faceValue(dVf_, facei)
                      + faceValue(dVfcorrectionValues, facei);

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

                    scalar dVff = faceValue(dVfcorrectionValues, facei);

                    dVff +=
                        sign(phi[fi])*fluidToPassThroughFace*sign(fluidToPassOn);

                    setFaceValue(dVfcorrectionValues, facei, dVff);

                    if (firstLoop)
                    {
                        checkIfOnProcPatch(facei);
                        correctedFaces.append(facei);
                    }
                }

                firstLoop = false;

                scalar alpha1New =
                (
                    alphaOld[celli]*rDeltaT  + Su[celli]
                  - netFlux(dVf_, celli)/Vi*rDeltaT
                  - netFlux(dVfcorrectionValues, celli)/Vi*rDeltaT
                )
                /
                (rDeltaT - Sp[celli]);

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


template<class SpType, class SuType>
void Foam::advection::isoAdvection::advect(const SpType& Sp, const SuType& Su)
{
    addProfilingInFunction(geometricVoF);
    DebugInFunction << endl;

    if (mesh_.topoChanging())
    {
        setProcessorPatches();
    }

    scalar advectionStartTime = mesh_.time().elapsedCpuTime();

    scalar rDeltaT = 1/mesh_.time().deltaTValue();

    // reconstruct the interface
    surf_->reconstruct();

    if(timeIndex_ < mesh_.time().timeIndex())
    {
        timeIndex_= mesh_.time().timeIndex();
        surf_->normal().oldTime() = surf_->normal();
        surf_->centre().oldTime() = surf_->centre();
    }

    // Initialising dVf with upwind values
    // i.e. phi[facei]*alpha1[upwindCell[facei]]*dt
    dVf_ = upwind<scalar>(mesh_, phi_).flux(alpha1_)*mesh_.time().deltaT();

    // Do the isoAdvection on surface cells
    timeIntegratedFlux();

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
      - fvc::surfaceIntegrate(dVf_)().primitiveField()*rDeltaT
    )/(rDeltaT - Sp.field());

    alpha1_.correctBoundaryConditions();

    // Adjust dVf for unbounded cells
    limitFluxes
    (
        Sp,
        Su
    );

    scalar maxAlphaMinus1 = gMax(alpha1In_) - 1;
    scalar minAlpha = gMin(alpha1In_);

    Info<< "isoAdvection: After conservative bounding: min(alpha) = "
        << minAlpha << ", max(alpha) = 1 + " << maxAlphaMinus1 << endl;

    // Apply non-conservative bounding mechanisms (clipping and snapping)
    // Note: We should be able to write out alpha before this is done!
    applyBruteForceBounding();

    advectionTime_ += (mesh_.time().elapsedCpuTime() - advectionStartTime);
    DebugInfo
        << "isoAdvection: time consumption = "
        << label(100*advectionTime_/(mesh_.time().elapsedCpuTime() + SMALL))
        << "%" << endl;

    alphaPhi_ = dVf_/mesh_.time().deltaT();
}


// ************************************************************************* //
