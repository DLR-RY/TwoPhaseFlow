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


#include "explicitInterfaceDiffFlux.H"
#include "zeroGradientFvPatchFields.H"

#include "mathematicalConstants.H"
#include "surfaceInterpolate.H"
#include "fvcDiv.H"
#include "fvcGrad.H"
#include "fvcSnGrad.H"

#include "plane.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::explicitInterfaceDiffFlux::explicitInterfaceDiffFlux
(
    const fvMesh& mesh
)
:
    mesh_(mesh)
{

}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::explicitInterfaceDiffFlux::diffusiveFlux
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

    const labelListList& stencil = exchangeFields.getStencil();

    Map<vector> mapCC =
        exchangeFields.getDatafromOtherProc(markedCells,mesh_.C());
    Map<scalar> mapPhi =
        exchangeFields.getDatafromOtherProc(markedCells,phi);

    forAll(markedCells,celli)
    {
        if (markedCells[celli] && mag(normal[celli]) != 0)
        {
            const scalar bcValue = boundaryValue[celli];
            vector n = normal[celli]/mag(normal[celli]);
            if (otherSide)
            {
                n *= -1;
            }

            scalar lowestAngle = -GREAT;
            scalar gradPhi = -GREAT;
            for (label i=1;i<stencil[celli].size();i++)
            {
                const label gblIdx = stencil[celli][i];
                const vector neiCC =
                    exchangeFields.getValue(mesh_.C(),mapCC,gblIdx);
                const vector dist = neiCC - centre[celli];
                scalar cosAngle = (dist/mag(dist)) & n;
                if (cosAngle > lowestAngle) // biggest possible value is 1
                {
                    lowestAngle = cosAngle;
                    scalar neiPhi =
                        exchangeFields.getValue(phi,mapPhi,gblIdx);
                    scalar deltaPhi = neiPhi-bcValue;
                    scalar d = mag(dist & n);
                    gradPhi = deltaPhi/d;
                }
            }

            if (gradPhi != -GREAT)
            {
                diffusiveFlux[celli] = gradPhi*gamma[celli];
            }

        }
    }

    return diffusiveFluxPtr;
}


Foam::tmp<Foam::fvScalarMatrix> Foam::explicitInterfaceDiffFlux::diffusiveFlux
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

    // is not necessarily update by reconstruct
    exchangeFields.setUpCommforZone(markedCells);

    const labelListList& stencil = exchangeFields.getStencil();

    Map<vector> mapCC =
        exchangeFields.getDatafromOtherProc(markedCells,mesh_.C());
    Map<scalar> mapPhi =
        exchangeFields.getDatafromOtherProc(markedCells,phi);


    forAll(markedCells,celli)
    {
        if (markedCells[celli] && mag(normal[celli]) != 0)
        {
            const scalar bcValue = boundaryValue[celli];
            const scalar area = mag(normal[celli]);
            vector n = normal[celli]/mag(normal[celli]);
            if (otherSide)
            {
                n *= -1;
            }

            scalar lowestAngle = -GREAT;
            scalar gradPhi = -GREAT;
            for (label i=1;i<stencil[celli].size();i++)
            {
                const label gblIdx = stencil[celli][i];
                const vector neiCC =
                    exchangeFields.getValue(mesh_.C(),mapCC,gblIdx);
                const vector dist = neiCC - centre[celli];
                scalar cosAngle = (dist/mag(dist)) & n;
                if (cosAngle > lowestAngle) // biggest possible value is 1
                {
                    lowestAngle = cosAngle;
                    scalar neiPhi =
                        exchangeFields.getValue(phi,mapPhi,gblIdx);
                    scalar deltaPhi = neiPhi-bcValue;
                    scalar d = mag(dist & n);
                    gradPhi = deltaPhi/d;
                }
            }

            if (gradPhi != -GREAT)
            {
                matrixSource[celli] -= gradPhi*gamma[celli]*area;
            }

        }
    }


    return diffusiveMatrixPtr;
}


void Foam::explicitInterfaceDiffFlux::diffusiveFlux
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

    // is not necessarily update by reconstruct
    exchangeFields.setUpCommforZone(markedCells);

    const labelListList& stencil = exchangeFields.getStencil();

    Map<vector> mapCC =
        exchangeFields.getDatafromOtherProc(markedCells,mesh_.C());
    Map<scalar> mapPhi =
         exchangeFields.getDatafromOtherProc(markedCells,phi);


    forAll(markedCells,celli)
    {
        if (markedCells[celli] && mag(normal[celli]) != 0)
        {
            const scalar bcValue = boundaryValue[celli];
            const scalar areaPerVol = mag(normal[celli])/mesh_.V()[celli];
            vector n = normal[celli]/mag(normal[celli]);
            if (otherSide)
            {
                n *= -1;
            }

            scalar lowestAngle = -GREAT;
            scalar gradPhi = -GREAT;
            for (label i=1;i<stencil[celli].size();i++)
            {
                const label gblIdx = stencil[celli][i];
                const vector neiCC =
                    exchangeFields.getValue(mesh_.C(),mapCC,gblIdx);
                const vector dist = neiCC - centre[celli];
                scalar cosAngle = (dist/mag(dist)) & n;
                if (cosAngle > lowestAngle) // biggest possible value is 1
                {
                    lowestAngle = cosAngle;
                    scalar neiPhi =
                        exchangeFields.getValue(phi,mapPhi,gblIdx);
                    scalar deltaPhi = neiPhi-bcValue;
                    scalar d = mag(dist & n);
                    gradPhi = deltaPhi/d;
                }
            }

            if (gradPhi != -GREAT)
            {
                matrixSource[celli] -= gradPhi*gamma[celli]*mag(normal[celli]);
                diffFlux[celli] += gradPhi*gamma[celli]*areaPerVol;
            }

        }
    }
}
