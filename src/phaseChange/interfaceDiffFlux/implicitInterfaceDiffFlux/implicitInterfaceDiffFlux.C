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


#include "implicitInterfaceDiffFlux.H"


Foam::label Foam::implicitInterfaceDiffFlux::getCellLabel
(
    const label& localIndex
)
{
    if (localIndex < mesh_.nCells())
    {
        return localIndex; // celli
    }
    else // cell next to it
    {
        const label faceI = localIndex
                            + mesh_.nInternalFaces()
                            - mesh_.nCells();

        const polyBoundaryMesh& pbm = mesh_.boundaryMesh();
        // Boundary face. Find out which face of which patch
        const label patchI = pbm.whichPatch(faceI);
        if (patchI < 0 || patchI >= pbm.size())
        {
            FatalErrorInFunction << "Cannot find patch for face " << faceI
                                << abort(FatalError);
        }
        const polyPatch& pp = pbm[patchI];
        const label patchFaceI = pp.whichFace(faceI);
        const label celli = pp.faceCells()[patchFaceI];
        return celli;
    }
}

void Foam::implicitInterfaceDiffFlux::getNeigbours
(
    const label celli,
    const vector faceCentre,
    const vector faceNormal,
    const volVectorField& centre,
    const volScalarField& phi,
    const Map<vector>& mapCC,
    const Map<scalar>& mapPhi,
    zoneDistribute& exchangeFields,
    DynamicField<scalar>& phiField,
    DynamicField<vector>& distField,
    DynamicField<label>& indicies
)
{
    const vector n = faceNormal/mag(faceNormal);
    const labelListList& stencil = exchangeFields.getStencil();
    for (label i=1;i<stencil[celli].size();++i)
    {
        const label gblIdx = stencil[celli][i];
        const vector neiCC = exchangeFields.getValue(centre,mapCC,gblIdx);
        const vector dist = neiCC - faceCentre;
        scalar cosAngle = (dist/mag(dist)) & n;
        if (cosAngle > 0.25) // roughly 75 deg
        {
            scalar neiPhi = exchangeFields.getValue(phi,mapPhi,gblIdx);
            phiField.append(neiPhi);
            distField.append(dist);
            indicies.append(gblIdx);
        }
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::implicitInterfaceDiffFlux::implicitInterfaceDiffFlux
(
    const fvMesh& mesh
)
:
    mesh_(mesh)
{

}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::implicitInterfaceDiffFlux::diffusiveFlux
(
    const boolList& markedCells,
    const volVectorField& normal,
    const volVectorField& centre,
    const volScalarField& phi,
    const volScalarField& gamma,
    const volScalarField& boundaryValue,
    const bool otherSide
)
{
    zoneDistribute& exchangeFields = zoneDistribute::New(mesh_);
    tmp<volScalarField> diffusiveFluxPtr
    (
        new volScalarField
        (
            IOobject
            (
                "diffusiveFlux",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar("init",phi.dimensions()*gamma.dimensions()/dimLength,0),
            calculatedFvPatchScalarField::typeName
        )
    );

    volScalarField& diffusiveFlux = diffusiveFluxPtr.ref();

    // is not necessarily update by reconstruct
    exchangeFields.setUpCommforZone(markedCells);


    Map<vector> mapCC =
        exchangeFields.getDatafromOtherProc(markedCells,mesh_.C());
    Map<scalar> mapPhi =
        exchangeFields.getDatafromOtherProc(markedCells,phi);

    DynamicField<scalar> neiPhiValues(100); // avoid resizing
    DynamicField<vector> neiDistValues(100); // avoid resizing
    DynamicField<label> neiIndicies(100); // avoid resizing

    forAll(markedCells,celli)
    {
        if (markedCells[celli] && mag(normal[celli]) != 0)
        {
            const scalar bcValue = boundaryValue[celli];
            scalar avgdPhidN = 0;
            scalar avgWeight = 0;
            neiPhiValues.clear();
            neiDistValues.clear();
            neiIndicies.clear();
            vector n = normal[celli]/mag(normal[celli]);
            if (otherSide)
            {
                n *= -1;
            }

            getNeigbours
            (
                celli,
                centre[celli],
                n,
                mesh_.C(),
                phi,
                mapCC,
                mapPhi,
                exchangeFields,
                neiPhiValues,
                neiDistValues,
                neiIndicies
            );

            forAll(neiPhiValues,i)
            {
                scalar phiNeiValue = neiPhiValues[i];
                vector dist = neiDistValues[i];
                scalar deltaPhi = phiNeiValue-bcValue;
                scalar d = mag(dist & n);
                scalar cosAngle = (dist/mag(dist)) & n;
                scalar weight = pow(mag(cosAngle),4);
                avgdPhidN += deltaPhi/d*weight;
                avgWeight += weight;
            }

            if (avgWeight != 0)
            {
                diffusiveFlux[celli] = avgdPhidN/avgWeight*gamma[celli];
            }

        }
    }

    return diffusiveFluxPtr;
}


Foam::tmp<Foam::fvScalarMatrix> Foam::implicitInterfaceDiffFlux::diffusiveFlux
(
    const boolList& markedCells,
    const volVectorField& normal,
    const volVectorField& centre,
    const volScalarField& phi,
    const volScalarField& gamma,
    const volScalarField& boundaryValue,
    const bool otherSide,
    const dimensionSet dim
)
{
    zoneDistribute& exchangeFields = zoneDistribute::New(mesh_);
    tmp<fvScalarMatrix> diffusiveMatrixPtr
    (
        new fvScalarMatrix(phi,dim)
    );

    fvScalarMatrix& diffusiveMatrix = diffusiveMatrixPtr.ref();
    scalarField& matrixSource = diffusiveMatrix.source();
    scalarField& matrixDiag = diffusiveMatrix.diag();

    // is not necessarily update by reconstruct
    exchangeFields.setUpCommforZone(markedCells);

    const globalIndex& gblNumbering = exchangeFields.globalNumbering();

    Map<vector> mapCC =
        exchangeFields.getDatafromOtherProc(markedCells,mesh_.C());
    Map<scalar> mapPhi =
        exchangeFields.getDatafromOtherProc(markedCells,phi);

    DynamicField<scalar> neiPhiValues(100); // avoid resizing
    DynamicField<vector> neiDistValues(100); // avoid resizing
    DynamicField<label> neiIndicies(100); // avoid resizing

    List<Map <FixedList<scalar,2> > > matrixDiagAndSource(Pstream::nProcs());


    forAll(markedCells,celli)
    {
        if (markedCells[celli] && mag(normal[celli]) != 0)
        {
            const scalar bcValue = boundaryValue[celli];
            const scalar area = mag(normal[celli]);
            scalar avgWeight = 0;
            neiPhiValues.clear();
            neiDistValues.clear();
            neiIndicies.clear();
            vector n = normal[celli]/mag(normal[celli]);
            if (otherSide)
            {
                n *= -1;
            }

            getNeigbours
            (
                celli,
                centre[celli],
                n,
                mesh_.C(),
                phi,
                mapCC,
                mapPhi,
                exchangeFields,
                neiPhiValues,
                neiDistValues,
                neiIndicies
            );

            forAll(neiPhiValues,i)
            {
                scalar cosAngle = (neiDistValues[i]/mag(neiDistValues[i])) & n;
                scalar weight = pow(mag(cosAngle),4);
                avgWeight += weight;
            }

            forAll(neiPhiValues,i)
            {
                vector dist = neiDistValues[i];
                scalar cosAngle = (dist/mag(dist)) & n;
                scalar weight = pow(mag(cosAngle),4);
                scalar weightDist = weight/mag(dist & n);
                scalar distFlux = weightDist/avgWeight*gamma[celli]*area;
                if (gblNumbering.isLocal(neiIndicies[i]))
                {
                    const label lclNeiCelli = getCellLabel
                    (
                        gblNumbering.toLocal(neiIndicies[i])
                    );
                    matrixDiag[lclNeiCelli] += distFlux;
                    matrixSource[lclNeiCelli] += distFlux*bcValue;
                }
                else
                {
                    const label gblIdx = neiIndicies[i];
                    const label procI = gblNumbering.whichProcID(gblIdx);
                    if (matrixDiagAndSource[procI].found(gblIdx))
                    {
                        FixedList<scalar,2>& matrixValues =
                             matrixDiagAndSource[procI][gblIdx];
                        matrixValues[0] += distFlux;
                        matrixValues[1] += distFlux*bcValue;
                    }
                    else
                    {
                        // other scalarField-entry is zero
                        FixedList<scalar,2> matrixValues(0.0);
                        matrixValues[0] = distFlux;
                        matrixValues[1] = distFlux*bcValue;
                        matrixDiagAndSource[procI].insert(gblIdx,matrixValues);
                    }
                }
            }
        }
    }

    if (Pstream::parRun())
    {
        Map<FixedList<scalar,2> >  MatrixCoeffFromOtherProc;

        PstreamBuffers pBufs(Pstream::commsTypes::nonBlocking);

        // Stream data into buffer
        for (label domain = 0; domain < Pstream::nProcs(); domain++)
        {
            if (domain != Pstream::myProcNo())
            {
                // Put data into send buffer
                UOPstream toDomain(domain, pBufs);

                toDomain << matrixDiagAndSource[domain];
            }
        }

        // wait until everything is written.
        pBufs.finishedSends();

        Map<FixedList<scalar,2>> neiMap;

        for (label domain = 0; domain < Pstream::nProcs(); domain++)
        {
            neiMap.clear();

            if (domain != Pstream::myProcNo())
            {
                // get data from send buffer
                UIPstream fromDomain(domain, pBufs);

                fromDomain >> neiMap;
                // append neiMap to boundaryAndProcData
                forAllConstIters(neiMap, iter)
                {
                    const label gblIdx = iter.key();

                    const label celli = getCellLabel
                    (
                        gblNumbering.toLocal(gblIdx)
                    );

                    matrixDiag[celli] += iter()[0];
                    matrixSource[celli] += iter()[1];


                }
            }
        }
    }

    return diffusiveMatrixPtr;
}


void Foam::implicitInterfaceDiffFlux::diffusiveFlux
(
    const boolList& markedCells,
    const volVectorField& normal,
    const volVectorField& centre,
    const volScalarField& phi,
    const volScalarField& gamma,
    const volScalarField& boundaryValue,
    const bool otherSide,
    volScalarField& diffFlux,
    fvScalarMatrix& phiSource
)
{
    zoneDistribute& exchangeFields = zoneDistribute::New(mesh_);

    scalarField& matrixSource = phiSource.source();
    scalarField& matrixDiag = phiSource.diag();

    // is not necessarily update by reconstruct
    exchangeFields.setUpCommforZone(markedCells);

    // const labelListList& stencil = exchangeFields.getStencil();
    const globalIndex& gblNumbering = exchangeFields.globalNumbering();

    Map<vector> mapCC =
        exchangeFields.getDatafromOtherProc(markedCells,mesh_.C());
    Map<scalar> mapPhi =
        exchangeFields.getDatafromOtherProc(markedCells,phi);

    DynamicField<scalar> neiPhiValues(100); // avoid resizing
    DynamicField<vector> neiDistValues(100); // avoid resizing
    DynamicField<label> neiIndicies(100); // avoid resizing

    List<Map <FixedList<scalar,2> > > matrixDiagAndSource(Pstream::nProcs());


    forAll(markedCells,celli)
    {
        if (markedCells[celli] && mag(normal[celli]) != 0)
        {
            const scalar bcValue = boundaryValue[celli];
            const scalar areaPerVol = mag(normal[celli])/mesh_.V()[celli];
            scalar avgdPhidN = 0;
            scalar avgWeight = 0;
            neiPhiValues.clear();
            neiDistValues.clear();
            neiIndicies.clear();
            vector n = normal[celli]/mag(normal[celli]);
            if (otherSide)
            {
                n *= -1;
            }

            getNeigbours
            (
                celli,
                centre[celli],
                n,
                mesh_.C(),
                phi,
                mapCC,
                mapPhi,
                exchangeFields,
                neiPhiValues,
                neiDistValues,
                neiIndicies
            );

            forAll(neiPhiValues,i)
            {
                scalar phiNeiValue = neiPhiValues[i];
                vector dist = neiDistValues[i];
                scalar deltaPhi = phiNeiValue-bcValue;
                scalar d = mag(dist & n);
                scalar cosAngle = (neiDistValues[i]/mag(neiDistValues[i])) & n;
                scalar weight = pow(mag(cosAngle),4);
                avgdPhidN += deltaPhi/d*weight;
                avgWeight += weight;
            }

            if (avgWeight != 0)
            {
                diffFlux[celli] = avgdPhidN/avgWeight*gamma[celli]*areaPerVol;
            }

            forAll(neiPhiValues,i)
            {
                vector dist = neiDistValues[i];
                scalar cosAngle = (dist/mag(dist)) & n;
                scalar weight = pow(mag(cosAngle),4);
                scalar weightDist = weight/mag(dist & n);
                scalar distFlux =
                    weightDist/avgWeight*gamma[celli]*mag(normal[celli]);
                if (gblNumbering.isLocal(neiIndicies[i]))
                {
                    const label lclNeiCelli = getCellLabel
                    (
                        gblNumbering.toLocal(neiIndicies[i])
                    );
                    matrixDiag[lclNeiCelli] += distFlux;
                    matrixSource[lclNeiCelli] += distFlux*bcValue;
                }
                else
                {
                    const label gblIdx = neiIndicies[i];
                    const label procI = gblNumbering.whichProcID(gblIdx);
                    if (matrixDiagAndSource[procI].found(gblIdx))
                    {
                        FixedList<scalar,2>& matrixValues =
                            matrixDiagAndSource[procI][gblIdx];
                        matrixValues[0] += distFlux;
                        matrixValues[1] += distFlux*bcValue;
                    }
                    else
                    {
                        // other scalarField-entry is zero
                        FixedList<scalar,2> matrixValues(0.0);
                        matrixValues[0] = distFlux;
                        matrixValues[1] = distFlux*bcValue;
                        matrixDiagAndSource[procI].insert(gblIdx,matrixValues);
                    }
                }
            }
        }
    }

    if (Pstream::parRun())
    {
        Map<FixedList<scalar,2> >  MatrixCoeffFromOtherProc;

        PstreamBuffers pBufs(Pstream::commsTypes::nonBlocking);

        // Stream data into buffer
        for (label domain = 0; domain < Pstream::nProcs(); domain++)
        {
            if (domain != Pstream::myProcNo())
            {
                // Put data into send buffer
                UOPstream toDomain(domain, pBufs);

                toDomain << matrixDiagAndSource[domain];
            }
        }

        // wait until everything is written.
        pBufs.finishedSends();

        Map<FixedList<scalar,2>> neiMap;

        for (label domain = 0; domain < Pstream::nProcs(); domain++)
        {
            neiMap.clear();

            if (domain != Pstream::myProcNo())
            {
                // get data from send buffer
                UIPstream fromDomain(domain, pBufs);

                fromDomain >> neiMap;
                // append neiMap to boundaryAndProcData
                forAllConstIters(neiMap, iter)
                {
                    const label gblIdx = iter.key();

                    const label celli = getCellLabel
                    (
                        gblNumbering.toLocal(gblIdx)
                    );

                    matrixDiag[celli] += iter()[0];
                    matrixSource[celli] += iter()[1];


                }
            }
        }
    }

}
// ************************************************************************* //
