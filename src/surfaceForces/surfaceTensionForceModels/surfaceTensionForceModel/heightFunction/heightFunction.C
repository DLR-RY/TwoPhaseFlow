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

#include "heightFunction.H"
#include "addToRunTimeSelectionTable.H"

#include "alphaContactAngleTwoPhaseFvPatchScalarField.H"
#include "mathematicalConstants.H"
#include "surfaceInterpolate.H"
#include "fvcDiv.H"
#include "fvcGrad.H"
#include "reconstructionSchemes.H"
#include "SortableList.H"
#include "leastSquareFitParabolid.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(heightFunction, 0);
    addToRunTimeSelectionTable(surfaceTensionForceModel,heightFunction, components);
}



void Foam::heightFunction::getStencilValues
(
    const Map<scalar>& MapAlpha,
    const labelList& stencil,
    DynamicField<scalar>& alphaValues
)
{
    alphaValues.clear();
    zoneDistribute& zDist = zoneDistribute::New(mesh_);
    forAll(stencil,i)
    {
        const label& gblIdx = stencil[i];
        alphaValues.append(zDist.getValue(alpha1_,MapAlpha,gblIdx));
    }
}

bool Foam::heightFunction::fullColumn(const scalar avgHeight,const scalar tol)
{
    return (avgHeight >= (1-tol) || avgHeight < tol);
}

void Foam::heightFunction::nextCell
(
    const label celli,
    const HFStencil::orientation orientation,
    twoDimFDStencil& HFCol
)
{
    label idx = stencilHF_.getCellLabel(celli,HFCol.nextCell(orientation));
    HFCol.status[orientation].gblIdx = idx;

}

void Foam::heightFunction::computeColumns
(
    const label dir,
    const Map<scalar>& mapAlpha,
    const label celli,
    const HFStencil::orientation orientation,
    twoDimFDStencil& HFCol
)
{
    const labelListList& stencil = stencilHF_.getStencil();
    const globalIndex& globalIdx = zoneDistribute::New(mesh_).globalNumbering();

    label& nextPos = HFCol.status[orientation].gblIdx;
    scalar columnHeightTol = 1e-6;

    nextCell
    (
        celli,
        orientation,
        HFCol
    );

    if (!globalIdx.isLocal(nextPos))
    {
        return;
    }

    scalar avgColVal = 0.5;
    DynamicField<scalar > alphaValues(100); // should be big enough avoids resizing

    for (label i = HFCol.status[orientation].iterI;i<7;i++) // move four times in both direction
    {
        if (globalIdx.isLocal(nextPos))
        {
            label localIdx = globalIdx.toLocal(nextPos);
            if (localIdx >= mesh_.nCells() || stencil[localIdx].size() == 0 )
            {
                // boundary face cannot compute height
                HFCol.status[orientation].iterI = i;
                HFCol.status[orientation].avgColVal = 0.5;
                return;
            }

            getStencilValues(mapAlpha,stencil[localIdx],alphaValues);

            HFCol.status[orientation].avgColVal = HFCol.addColumnHeight(alphaValues);

            if (!fullColumn(HFCol.status[orientation].avgColVal ,columnHeightTol))
            {
                HFCol.status[orientation].iterI = i;

                nextCell
                (
                    localIdx,
                    orientation,
                    HFCol
                );
            }
        }

        if (fullColumn(HFCol.status[orientation].avgColVal,columnHeightTol) )
        {
            return;
        }

    }
}

void Foam::heightFunction::exchangeStencils
(
    const List<List<twoDimFDStencil>>& sendStencil,
    List<twoDimFDStencil>& recvStencil
)
{
    PstreamBuffers pBufs(Pstream::commsTypes::nonBlocking);

    // Stream data into buffer
    for (label domain = 0; domain < Pstream::nProcs(); domain++)
    {
        if (domain != Pstream::myProcNo())
        {
            // Put data into send buffer
            UOPstream toDomain(domain, pBufs);

            toDomain << sendStencil[domain];
        }
    }

    // wait until everything is written.
    pBufs.finishedSends();

    for (label domain = 0; domain < Pstream::nProcs(); domain++)
    {
        if (domain != Pstream::myProcNo())
        {
            List<twoDimFDStencil> stencilFromDomain;
            // get data from send buffer
            UIPstream fromDomain(domain, pBufs);

            fromDomain >> stencilFromDomain;
            recvStencil.append(stencilFromDomain);

        }
    }
}

Foam::scalar Foam::heightFunction::calcCurvature
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
        return (2*fit[2] + 2*fit[3] + 2*fit[2]*sqr(fit[1]) + 2*fit[3]*sqr(fit[0]) - 2*fit[4]*fit[0]*fit[1])/pow(1+sqr(fit[0])+sqr(fit[1]),1.5);
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::heightFunction::heightFunction
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
    ),
    mesh_(alpha1.mesh()),
    stencilHF_(alpha1.mesh(),dict.lookupOrDefault<scalar>("angleTol",0.001)),
    IFRegion_(alpha1.mesh()),
    twoDim_(false)

{
    label dimensions = 0;
    for (direction cmpt = 0; cmpt < vector::nComponents; cmpt++)
    {
        dimensions += pos0(mesh_.geometricD()[cmpt]) * mesh_.geometricD()[cmpt];
    }

    twoDim_ = (dimensions == 2);

    if (dimensions <= 1)
    {
        FatalErrorInFunction
        << "The Height Function Method only works on grids with atleast " << nl
        << "two dimensions" << nl
        << endl
        << exit(FatalError);
    }
}

// * * * * * * * * * * * * * * Public Access Member Functions  * * * * * * * * * * * * * * //

void Foam::heightFunction::correctContactAngle
(
    surfaceVectorField::Boundary& nHatb,
    surfaceVectorField::Boundary& heightFunctionf
)
{
    scalar convertToRad = Foam::constant::mathematical::pi/180.0;

    const fvMesh& mesh = alpha1_.mesh();
    const volScalarField::Boundary& abf = alpha1_.boundaryField();

    const fvBoundaryMesh& boundary = mesh.boundary();

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

            fvsPatchVectorField& nHatp = nHatb[patchi];
            const scalarField theta
            (
                convertToRad*acap.theta(U_.boundaryField()[patchi], nHatp)
            );

            const vectorField nf
            (
                boundary[patchi].nf()
            );

            // Reset nHatp to correspond to the contact angle

            const scalarField a12(nHatp & nf);
            const scalarField b1(cos(theta));

            scalarField b2(nHatp.size());
            forAll(b2, facei)
            {
                b2[facei] = cos(acos(a12[facei]) - theta[facei]);
            }

            const scalarField det(1.0 - a12*a12);

            scalarField a((b1 - a12*b2)/det);
            scalarField b((b2 - a12*b1)/det);

            nHatp = a*nf + b*nHatp;
            nHatp /= (mag(nHatp) + deltaN_.value());

            acap.gradient() = (nf & nHatp)*mag(heightFunctionf[patchi]);
            acap.evaluate();
        }
    }
}


void Foam::heightFunction::correct()
{
    deltaFunctionModel_->correct();

    const fvMesh& mesh = alpha1_.mesh();
    const surfaceVectorField& Sf = mesh.Sf();

    reconstructionSchemes& surf =
        mesh.lookupObjectRef<reconstructionSchemes>("reconstructionScheme");

    // can also be an isosurface
    surf.reconstruct(false);

    const boolList& interfaceCells = surf.interfaceCell();
    const volVectorField& faceCentre = surf.centre();
    const volVectorField& faceNormal = surf.normal();
    const globalIndex& globalIdx = zoneDistribute::New(mesh_).globalNumbering();

    boolList nextToInterface(mesh.nCells(),false);

    scalar deltaX = mag(mesh_.delta())().average().value();

    volScalarField cellDistField
    (
        IOobject
        (
            "cellDistField",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("scalar", dimless, 0),
        "calculated"   //alpha1.boundaryField().types()
    );

    IFRegion_.markCellsNearSurf(interfaceCells,5,nextToInterface,cellDistField);

    if (mesh_.time().outputTime())
    {
        cellDistField.write();
    }
    // updateStencil;
    stencilHF_.updateStencil(nextToInterface);

    const boolList& isCuboid = stencilHF_.isCuboid();

    volScalarField isCuboidField
    (
        IOobject
        (
            "isCuboidField",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("scalar", dimless, 0),
        "calculated"   //alpha1.boundaryField().types()
    );

    forAll(isCuboid,i)
    {
        if (isCuboid[i])
        {
            isCuboidField[i] = 1;
        }
        else
        {
            isCuboidField[i] = 0;
        }
    }

    if (mesh_.time().outputTime())
    {
        isCuboidField.write();
    }

    Map<scalar> MapAlpha = stencilHF_.getDatafromOtherProc
    (
        nextToInterface,
        alpha1_
    );

    volScalarField foundHeightField
    (
        IOobject
        (
            "foundHeightField",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("scalar", dimless, 0),
        "calculated"   //alpha1.boundaryField().types()
    );
    label test = twoDim_ ? 2 : 3;
    Vector<label> geomDir = mesh_.geometricD();
    label dirs = twoDim_ ? 2 : 3;

    twoDimStencilMap parallelStencil;
    List<List<twoDimFDStencil>> sendStencil(Pstream::nProcs());

    DynamicField<scalar > alphaValues(100); // should be big enough avoids resizing

    const labelListList& stencil = stencilHF_.getStencil();

    forAll(interfaceCells,celli)
    {
        vector n = faceNormal[celli];
        if (interfaceCells[celli] && mag(n) != 0)
        {
            if (isCuboid[celli])
            {
                n /= mag(n);

                SortableList<scalar> sortedDirs(0);
                for (int dirI=0;dirI<3;dirI++)
                {
                    if (geomDir[dirI] != -1)
                        sortedDirs.append(mag(n[dirI]));
                }
                sortedDirs.reverseSort();
                foundHeightField[celli] = 0;

                bool success = false;
                for (int dirI=0;dirI<dirs;dirI++)
                {
                    label direction = sortedDirs.indices()[dirI];
                    twoDimFDStencil cols(twoDim_,direction,globalIdx.toGlobal(celli));
                    cols.status.first().iterI = 1;
                    cols.status.second().iterI = 1;

                    getStencilValues(MapAlpha,stencil[celli],alphaValues);
                    cols.addColumnHeight(alphaValues);

                    computeColumns
                    (
                        direction,
                        MapAlpha,
                        celli,
                        HFStencil::orientation::pos,
                        cols
                    );

                    computeColumns
                    (
                        direction,
                        MapAlpha,
                        celli,
                        HFStencil::orientation::neg,
                        cols
                    );

                    if (cols.foundHeight())
                    {
                        foundHeightField[celli] = 1;
                        K_[celli] = cols.calcCurvature(deltaX);
                        success = true;
                        break;
                    }

                    if (!globalIdx.isLocal(cols.status[0].gblIdx))
                    {
                        parallelStencil.insert
                        (
                            Vector2D<label>(celli,direction),
                            cols
                        );
                        label procI = globalIdx.whichProcID(cols.status[0].gblIdx);
                        sendStencil[procI].append(cols);
                    }
                    if (!globalIdx.isLocal(cols.status[1].gblIdx))
                    {
                        parallelStencil.insert
                        (
                            Vector2D<label>(celli,direction),
                            cols
                        );
                        label procI = globalIdx.whichProcID(cols.status[1].gblIdx);
                        sendStencil[procI].append(cols);
                    }

                }
            }
            else
            {
                // fitParaboloid
            }
        }
        else
        {
            // fitParaboloid
            K_[celli] = 0;
        }

    }

    List<twoDimFDStencil> recvStencil;
    exchangeStencils(sendStencil,recvStencil);

    sendStencil.clear();
    sendStencil.resize(Pstream::nProcs());

    forAll(recvStencil,i)
    {
        twoDimFDStencil& HFcols = recvStencil[i];
        HFcols.heights() = 0.0;

        if (globalIdx.isLocal(HFcols.status[0].gblIdx))
        {
            const label localIdx = globalIdx.toLocal(HFcols.status[0].gblIdx);
            HFcols.status[1].gblIdx = -1;

            computeColumns
            (
                HFcols.direction(),
                MapAlpha,
                localIdx,
                HFStencil::orientation::pos,
                HFcols
            );

            if
            (
                HFcols.status[0].avgColVal < 1e-6
             || HFcols.status[0].avgColVal >  1-1e-6
            )
            {
                label origProc = globalIdx.whichProcID(HFcols.gblcelli());
                sendStencil[origProc].append(HFcols);
            }
        }

        if (globalIdx.isLocal(HFcols.status[1].gblIdx))
        {
            const label localIdx = globalIdx.toLocal(HFcols.status[1].gblIdx);
            HFcols.status[0].gblIdx = -1;

            computeColumns
            (
                HFcols.direction(),
                MapAlpha,
                localIdx,
                HFStencil::orientation::neg,
                HFcols
            );

            if
            (
                HFcols.status[0].avgColVal < 1e-6
             || HFcols.status[0].avgColVal >  1-1e-6
            )
            {
                label origProc = globalIdx.whichProcID(HFcols.gblcelli());
                sendStencil[origProc].append(HFcols);
            }
        }

        // send back otherProc
    }

    List<twoDimFDStencil> sendBackToOriginalProcess;
    exchangeStencils(sendStencil,sendBackToOriginalProcess);

    forAll(sendBackToOriginalProcess,i)
    {
        twoDimFDStencil& neiCols = sendBackToOriginalProcess[i];
        if (neiCols.status[0].gblIdx != -1)
        {
            label celli = globalIdx.toLocal(neiCols.gblcelli());
            Vector2D<label> celliAndDir(celli,neiCols.direction());
            if (parallelStencil.found(celliAndDir))
            {
                parallelStencil[celliAndDir].status[0] = neiCols.status[0];
                parallelStencil[celliAndDir].heights() += neiCols.heights();
            }
        }
        if (neiCols.status[1].gblIdx != -1)
        {
            label celli = globalIdx.toLocal(neiCols.gblcelli());
            Vector2D<label> celliAndDir(celli,neiCols.direction());
            if (parallelStencil.found(celliAndDir))
            {
                parallelStencil[celliAndDir].status[1] = neiCols.status[1];
                parallelStencil[celliAndDir].heights() += neiCols.heights();
            }
        }
    }

    forAllIters(parallelStencil,iter)
    {
        if (iter().foundHeight())
        {
            label celli = globalIdx.toLocal(iter().gblcelli());
            foundHeightField[celli] = 1;
            K_[celli] = iter().calcCurvature(deltaX);
        }

    }

    // fitParaboloid

    Vector<label> explicitDim(1,1,-1);

    label nDims = 0;

    forAll(mesh.geometricD(),i)
    {
        if (mesh.geometricD()[i] == 1)
        {
            nDims++;
        }
    }

    if (nDims == 2)
    {
        explicitDim.y() = -1;
    }

    leastSquareFitParabolid paraboloid(geomDir,explicitDim);

    zoneDistribute& zDist = zoneDistribute::New(mesh_);

    zDist.setUpCommforZone(interfaceCells);

    Map<Field <vector > > mapCentres = zDist.getFields(interfaceCells,faceCentre);
    Map<Field <vector > > mapNormal = zDist.getFields(interfaceCells,faceNormal);

    forAll(interfaceCells,cellI)
    {
        if (interfaceCells[cellI] && foundHeightField[cellI] == 0)
        {
            if (mag(faceNormal[cellI]) == 0)
            {
                K_[cellI] = 0;
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
                K_[cellI] =  calcCurvature(paraboloid.fitParaboloid(c,n,centres,weight));
            }
            else
            {
                K_[cellI] = 0;
            }
        }
    }

    if (mesh_.time().outputTime())
    {
        foundHeightField.write();
    }
    nextToInterface = false;

    IFRegion_.markCellsNearSurf(interfaceCells,1,nextToInterface,cellDistField);

    zDist.setUpCommforZone(nextToInterface);

    Map<Field<scalar>> mapCurv = zDist.getFields(nextToInterface,K_);
    mapCentres = zDist.getFields(nextToInterface,faceCentre);

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
                if (centreField[i] != vector::zero)
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

    Kf_ = fvc::interpolate(K_);

}



// ************************************************************************* //
