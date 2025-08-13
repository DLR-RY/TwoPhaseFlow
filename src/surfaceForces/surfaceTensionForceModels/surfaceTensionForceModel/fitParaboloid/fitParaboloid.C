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

#include "fitParaboloid.H"
#include "addToRunTimeSelectionTable.H"

#include "alphaContactAngleTwoPhaseFvPatchScalarField.H"
#include "mathematicalConstants.H"
#include "surfaceInterpolate.H"
#include "fvcDiv.H"
#include "fvcGrad.H"
#include "fvcSnGrad.H"
#include "fvc.H"

#include "plane.H"
#include "interpolationCellPoint.H"
#include "tensor2D.H"

#include "reconstructionSchemes.H"
#include "leastSquareFitParabolid.H"
#include "cutFacePLIC.H"
#include "cutCellIso.H"
#include "reconstructedDistanceFunction.H"
#include "processorPolyPatch.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(fitParaboloid, 0);
    addToRunTimeSelectionTable(surfaceTensionForceModel,fitParaboloid, components);
}


Foam::scalar Foam::fitParaboloid::calcCurvature
(
    const scalarField& fit
)
{
    if (fit.size() == 2)
    {
        return (2*fit[1])/pow(1+sqr(fit[0]),1.5);
    }
    else
    {
        // c0*x + c1*y  + c2*x^2 + c3*y^2  + c4*x*y
        // curvature at x=0, y=0
        // 1st derivative = c1, c2
        // 2nd derivative = 2*c3, 2*c4, c5
        // equation: Jibben et al, A Paraboloid Fitting Technique for
        // Calculating Curvature from Piecewise-Linear
        // Interface Reconstructions on 3D
        // Unstructured Meshes
        return (2*fit[2] + 2*fit[3] + 2*fit[2]*sqr(fit[1])
                + 2*fit[3]*sqr(fit[0]) - 2*fit[4]*fit[0]*fit[1])
                /pow(1+sqr(fit[0])+sqr(fit[1]),1.5);
    }
}

Foam::vectorField Foam::fitParaboloid::getFaceCentres
(
    const globalIndex& globalNumbering,
    const DynamicList < label >& stencil,
    const volVectorField& faceCentres,
    const Map < vector >& map
)
{
    const label nCells = faceCentres.mesh().nCells();

    DynamicList<vector> centres(stencil.size());
    forAll(stencil, i)
    {
        const label gIdx = stencil[i];
        if (globalNumbering.isLocal(gIdx))
        {

            // also includes boundary conditions
            const label localCellI = globalNumbering.toLocal( gIdx);
            if (localCellI < nCells)
            {
                if (mag(faceCentres[localCellI]) != 0)
                {
                    centres.append(faceCentres[localCellI]);
                }
            }
            else
            {
                // face is saved in map
                centres.append(map[gIdx]);
            }
        }
        else
        {
            // face is saved in map
            centres.append(map[gIdx]);

        }
    }

    return vectorField(centres);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fitParaboloid::fitParaboloid
(
    const dictionary& dict,
    const volScalarField& alpha1,
    const surfaceScalarField& phi,
    const volVectorField& U
)
:
    surfaceTensionForceModel
    (
        typeName,
        dict,
        alpha1,
        phi,
        U
    ),
    deltaN_
    (
        "deltaN",
        1e-8/pow(average(alpha1.mesh().V()), 1.0/3.0)
    )
{

}


// * * * * * * * * * * * * * * Public Access Member Functions  * * * * * * * * * * * * * * //

void Foam::fitParaboloid::correctContactAngle
(
    surfaceVectorField::Boundary& nHatb,
    surfaceVectorField::Boundary& gradAlphaf
)
{

}


void Foam::fitParaboloid::correctContactAngle
(
    volVectorField& normal,
    volVectorField& centre
)
{
    scalar convertToRad = Foam::constant::mathematical::pi/180.0;

    const fvMesh& mesh = alpha1_.mesh();

    // check if face is cut
    cutFacePLIC cutFace(mesh);

    const volScalarField::Boundary& abf = alpha1_.boundaryField();
    volVectorField::Boundary& cbf = centre.boundaryFieldRef();
    volVectorField::Boundary& nbf = normal.boundaryFieldRef();

    const fvBoundaryMesh& boundary = mesh.boundary();

    // we need a surfaceVectorField to compute theta
    surfaceVectorField normalf(fvc::interpolate(normal));

    forAll(boundary, patchi)
    {
        if (isA<alphaContactAngleTwoPhaseFvPatchScalarField>(abf[patchi]))
        {
            forAll(normalf.boundaryFieldRef()[patchi],i)
            {
                const label celli = boundary[patchi].faceCells()[i];
                vector n = normal[celli];
                if (mag(n) != 0)
                {
                    n /= mag(n);
                    normalf.boundaryFieldRef()[patchi][i] = n;
                }
            }
        }
    }

    forAll(boundary, patchi)
    {
        if (isA<alphaContactAngleTwoPhaseFvPatchScalarField>(abf[patchi]))
        {
            alphaContactAngleTwoPhaseFvPatchScalarField& acap =
                const_cast<alphaContactAngleTwoPhaseFvPatchScalarField&>
                (
                    refCast<const alphaContactAngleTwoPhaseFvPatchScalarField>
                    (
                        abf[patchi]
                    )
                );

            fvsPatchVectorField& nHatp = normalf.boundaryFieldRef()[patchi];
            const scalarField theta
            (
                convertToRad*acap.theta(U_.boundaryField()[patchi], nHatp)
            );

            const vectorField nf
            (
                boundary[patchi].nf()
            );

            // Reset nHatp to correspond to the contact angle
            forAll(nbf[patchi],i)
            {
                const label celli = boundary[patchi].faceCells()[i];
                const label faceI = boundary[patchi].start() + i;
                vector n = normal[celli];
                if (mag(n) != 0)
                {
                    n /= mag(n);
                    label cutStatus = cutFace.calcSubFace
                    (
                        faceI,
                        n,
                        centre[celli]
                    );

                    if (cutStatus == 0)
                    {
                        const point cutEdgeCentre =
                             average(cutFace.surfacePoints());

                        // project Normal on the face
                        vector projN = (tensor::I - nf[i]*nf[i]) & -n;

                        // normalise
                        projN /= mag(projN) +deltaN_.value();

                        vector nTheta =
                            sin(theta[i])*nf[i] + cos(theta[i])*projN;

                        scalar proJDist =
                            mag((boundary[patchi].Cf()[i] - centre[celli]) & nf[i]);

                        // should point outside of the domain
                        cbf[patchi][i] =
                            cutEdgeCentre
                          + nTheta/boundary[patchi].deltaCoeffs()[i];

                        nbf[patchi][i] = normal[celli];
                    }
                }
                else
                {
                    cbf[patchi][i] = vector::zero;
                }
            }
            acap.evaluate();
        }
    }
}


void Foam::fitParaboloid::correct()
{
    deltaFunctionModel_->correct();

    const fvMesh& mesh = alpha1_.mesh();

    reconstructionSchemes& surf = mesh.lookupObjectRef<reconstructionSchemes>("reconstructionScheme");
    // can also be an isosurface
    surf.reconstruct(false);

    zoneDistribute& exchangeFields = zoneDistribute::New(mesh);

    boolList interfaceCells = surf.interfaceCell();

    volVectorField faceCentre = surf.centre();

    volVectorField faceNormal = surf.normal();

    Vector<label> geomDir = mesh.geometricD();
    Vector<label> explicitDim(1,1,-1);

    label nDims = 0;

    forAll(geomDir,i)
    {
        if (geomDir[i] == 1)
        {
            nDims++;
        }
    }

    if (nDims == 2)
    {
        explicitDim.y() = -1;
    }

    leastSquareFitParabolid paraboloid(geomDir,explicitDim);

    Map<Field <vector > > mapCentres;
    Map<Field <vector > > mapNormal;

    exchangeFields.setUpCommforZone(interfaceCells);

    mapCentres = exchangeFields.getFields(interfaceCells,faceCentre);
    mapNormal = exchangeFields.getFields(interfaceCells,faceNormal);

    boolList nextToInterface(mesh.nCells(),false);
    const globalIndex globalNumbering = exchangeFields.globalNumbering();

    forAll(interfaceCells,cellI)
    {
        if (interfaceCells[cellI])
        {
            if (mag(faceNormal[cellI]) == 0)
            {
                K_[cellI] = 0;
                interfaceCells[cellI] = false;
                nextToInterface[cellI] = true;
                continue;
            }
            vector n = faceNormal[cellI]/mag(faceNormal[cellI]);
            point c = faceCentre[cellI];

            const vectorField& neiNormal =  mapNormal[cellI];
            const vectorField& neiCentre =  mapCentres[cellI];

            DynamicField< vector > centres(neiCentre[cellI].size());
            DynamicField< scalar > weight(neiNormal[cellI].size());

            forAll(neiNormal,i)
            {
                if (mag(neiNormal[i]) != 0)
                {
                    centres.append(neiCentre[i]);
                    weight.append(pow(mag(neiNormal[i]),0.25));
                }
            }

            if (centres.size() >= paraboloid.nCoeffs())
            {
                K_[cellI] =
                    calcCurvature(paraboloid.fitParaboloid(c,n,centres,weight));
            }
            else
            {
                K_[cellI] = 0;
            }

            forAll(exchangeFields.getStencil()[cellI],i)
            {
                const label gblIdx = exchangeFields.getStencil()[cellI][i];
                if (globalNumbering.isLocal(gblIdx))
                {
                    const label idx = globalNumbering.toLocal(gblIdx);
                    if (idx < mesh.nCells())
                    {
                        nextToInterface[idx] = true;
                    }
                }

            }

        }
        else
        {
            K_[cellI] = 0;
        }
    }

    K_.correctBoundaryConditions();

    Map<Field <scalar > > mapCurv;

    exchangeFields.setUpCommforZone(nextToInterface);

    mapCurv = exchangeFields.getFields(nextToInterface,K_);
    mapCentres = exchangeFields.getFields(nextToInterface,faceCentre);

    forAll(nextToInterface,celli)
    {
        if (nextToInterface[celli] && !interfaceCells[celli])
        {
            const point cc = mesh.C()[celli];
            scalar smallDist = GREAT;
            label smallestDistIdx = -1;
            Field < vector> centreField = mapCentres[celli];
            forAll(centreField,i)
            {
                if (centreField[i] != vector::zero && i > 1)
                {
                    scalar dist = mag(cc-centreField[i]);
                    if (smallDist > dist)
                    {
                        smallDist = dist;
                        smallestDistIdx = i;
                    }
                }
            }
            K_[celli] = mapCurv[celli][smallestDistIdx];
        }
    }

    const labelList& owner = mesh.owner();
    const labelList& neighbour = mesh.neighbour();

    forAll(Kf_,faceI)
    {
        const label own = owner[faceI];
        const label nei = neighbour[faceI];

        // both are true
        if (interfaceCells[own] && interfaceCells[nei])
        {
            scalar weight = 0.5;
            Kf_[faceI] = weight*K_[own] + (1-weight)*K_[nei];
        }
        else if (interfaceCells[own] != interfaceCells[nei]) // one is true
        {
            if (interfaceCells[own])
            {
                Kf_[faceI] = K_[own];
            }
            else
            {
                Kf_[faceI] = K_[nei];
            }
        }
        else // both are false
        {
            Kf_[faceI] = 0;
        }
    }

    const polyBoundaryMesh& patches = mesh.boundaryMesh();


    PstreamBuffers pBufs(Pstream::commsTypes::nonBlocking);

    List<DynamicList<label >>  procPatchLabels(patches.size());
    List<DynamicList<scalar >> procPatchValue(patches.size());

    forAll(Kf_.boundaryFieldRef(),patchi)
    {
        const polyPatch& patch = patches[patchi];
        bool isprocPatch = isA<processorPolyPatch>(patches[patchi]);
        // you can not write on empty patches
        if (!isA<emptyPolyPatch>(patches[patchi]))
        {
            forAll(patch.faceCells(),i)
            {
                const label celli = patch.faceCells()[i];
                if (interfaceCells[celli])
                {
                    Kf_.boundaryFieldRef()[patchi][i] = K_[celli];
                    if (isprocPatch)
                    {
                        procPatchLabels[patchi].append(i);
                        procPatchValue[patchi].append(K_[celli]);
                    }
                }
                else
                {
                    Kf_.boundaryFieldRef()[patchi][i] = 0;
                }
            }
        }
    }


    if (Pstream::parRun())
    {
        PstreamBuffers pBufs(Pstream::commsTypes::nonBlocking);

        // Send
        forAll(procPatchLabels, patchi)
        {
            if (isA<processorPolyPatch>(patches[patchi]))
            {
            const processorPolyPatch& procPatch =
                refCast<const processorPolyPatch>(patches[patchi]);

            UOPstream toNbr(procPatch.neighbProcNo(), pBufs);

            toNbr << procPatchLabels[patchi] << procPatchValue[patchi];
            }

        }

        pBufs.finishedSends();

        // Receive and combine
        forAll(procPatchLabels, patchi)
        {
            if (isA<processorPolyPatch>(patches[patchi]))
            {

                const processorPolyPatch& procPatch =
                    refCast<const processorPolyPatch>(patches[patchi]);

                UIPstream fromNeighb(procPatch.neighbProcNo(), pBufs);
                List<label> faceIDs;
                List<scalar> nbrdCurvs;

                fromNeighb >> faceIDs >> nbrdCurvs;
                scalarField& localCurv = Kf_.boundaryFieldRef()[patchi];

                forAll(faceIDs, i)
                {
                    const label facei = faceIDs[i];
                    if (localCurv[facei] !=  0)
                    {
                        localCurv[facei] =
                            0.5*(localCurv[facei] + nbrdCurvs[i]);
                    }
                    else
                    {
                        localCurv[facei] =  nbrdCurvs[i];
                    }
                }

            }
        }
    }
}



// ************************************************************************* //
