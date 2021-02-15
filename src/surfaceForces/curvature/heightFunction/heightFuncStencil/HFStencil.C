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

#include "HFStencil.H"
#include "dummyTransform.H"
#include "emptyPolyPatch.H"
#include "processorPolyPatch.H"
#include "syncTools.H"
#include "wedgePolyPatch.H"
#include "unitConversion.H"


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::HFStencil::findCuboids()
{
    const cellList& cells = mesh_.cells();
    const vectorField& faceCentre = mesh_.faceCentres();

    scalar angleError0Deg = cos(degToRad(toleranceAngle_));
    scalar angleError90Deg = cos(degToRad(90-toleranceAngle_));

    forAll(cells,celli)
    {
        const cell& c = cells[celli];
        if (c.nFaces() != 6)
        {
            isCuboid_[celli] = false;
            continue;
        }

        const vector cc = mesh_.C()[celli];
        bool cuboid = true;

        // find biggest x direction
        label maxDir = -1;
        scalar maxValX = -10;
        vector n(0,0,0);
        forAll(c,i)
        {
            const label& faceI = c[i];
            vector dist = faceCentre[faceI]-cc;
            dist /= mag(dist);
            if (dist.x() >= maxValX)
            {
                maxDir = i;
                maxValX = dist.x();
                n = dist;
            }
        }

        forAll(c,i)
        {
            const label& faceI = c[i];
            vector dist = faceCentre[faceI]-cc;
            dist /= mag(dist);
            scalar angle = mag(n & dist);
            // angle between 0 and zero
            if (angle >=  angleError90Deg && angle <= angleError0Deg)
            {
                cuboid = false;
                break;
            }

        }
        if (cuboid)
        {
            isCuboid_[celli] = true;
        }
    }

    const polyBoundaryMesh& pbm = mesh_.boundaryMesh();
    forAll(pbm, patchi) // mark on point patches
    {
        const polyPatch& pp = pbm[patchi];

        // exclude empty and wedgepatches
        if
        (
            !isA<processorPolyPatch>(pp)
            && !isA<wedgePolyPatch>(pp)
            && !isA<emptyPolyPatch>(pp)
        )
        {
            forAll(pp.faceCells(),i)
            {
                const label pCelli = pp.faceCells()[i];
                isCuboid_[pCelli] = false;
            }
        }
    }
}


Foam::label Foam::HFStencil::categorizePoint
(
    const scalar pos
)
{

    label relativePos = 1; // in the middle

    if (pos <= -(sin(degToRad(toleranceAngle_))))
    {
        relativePos = 0; // in the left
    }
    else if (pos >= sin(degToRad(toleranceAngle_)))
    {
        relativePos = 2; // in the right
    }

    return relativePos;

}


Foam::label Foam::HFStencil::calcPosInStencil
(
    const vector& cc,
    const vector& neic
)
{
    // check if two Dim
    Vector< label > structAddress(0,0,0);
    //i = j = k = 0;
    // i is x
    // j is y
    // k is z
    // if (twoDim_ == -1)
    // {
        // 3D case
    vector distCellToCell = neic - cc;
    if (mag(distCellToCell))
        distCellToCell /= mag(distCellToCell);
    // Info << "distCellToCell " << distCellToCell << endl;

    structAddress.x() = categorizePoint(distCellToCell.x());
    structAddress.y() = categorizePoint(distCellToCell.y());
    structAddress.z() = categorizePoint(distCellToCell.z());

    // Info << "structAddress " << structAddress << endl;

    // set k to Zero only use i, j
    if (twoDim_)
    {
        bool foundFirstDir = false;
        for (label geomDir=0;geomDir<3;geomDir++)
        {
            if (mesh_.geometricD()[geomDir] == 1)
            {
                if (!foundFirstDir)
                {
                    structAddress[0] = structAddress[geomDir];
                    foundFirstDir = true;
                }
                else
                {
                    structAddress[1] = structAddress[geomDir];
                }

            }
        }
        structAddress[2] = 0;
    }

    return structAddress.x() + 3*structAddress.y() + structAddress.z()*9;
}


void Foam::HFStencil::sortStencilToIJKFormat
(
    const label celli,
    const labelList& CPCstencil,
    const vectorField& cellCentres,
    const vector cc
)
{
    // sort
    if (twoDim_)
    {
        stencilHF_[celli].setSize(9);
    }
    else
    {
        stencilHF_[celli].setSize(27);
    }

    stencilHF_[celli] = -1;

    forAll(CPCstencil,i)
    {
        label labelInIJK = calcPosInStencil(mesh_.C()[celli],cellCentres[i]);
        stencilHF_[celli][labelInIJK] = CPCstencil[i];
    }
}


Foam::label Foam::HFStencil::arrayToList3D(const label i,const label j,const label k)
{
    // assumes that i,j and k are between 0 and 2
    return (i + 3*j + 9*k);
}


Foam::label Foam::HFStencil::getCellLabel
(
    const label celli,
    const Vector<label>& posInStencil
)
{
    if (isCuboid_[celli])
    {
        if (posInStencil.z() != 0 && stencilHF_[celli].size() == 9)
        {
            return -1;
        }
        label idx = arrayToList3D
        (
            posInStencil.x(),
            posInStencil.y(),
            posInStencil.z()
        );
        return stencilHF_[celli][idx];
    }
    else
    {
        return -1;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::HFStencil::HFStencil(const fvMesh& mesh,const scalar angleTolerance)
:
    mesh_(mesh),
    stencilHF_(mesh.nCells()),
    isCuboid_(mesh.nCells(), false),
    twoDim_(false),
    toleranceAngle_(angleTolerance)
{
    label dimensions = 0;
    for (Foam::direction cmpt = 0; cmpt < vector::nComponents; cmpt++)
    {
        dimensions += pos0(mesh_.geometricD()[cmpt]) * mesh_.geometricD()[cmpt]; // why the last one
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

    findCuboids();

}


void Foam::HFStencil::updateStencil(const boolList& nextToInterface)
{
    if (mesh_.topoChanging())
    {
        // Introduced resizing to cope with changing meshes
        isCuboid_.resize(mesh_.nCells());
        isCuboid_ = false;
        findCuboids();
        stencilHF_ = labelListList(mesh_.nCells());
    }

    zoneDistribute& exchangeFields = zoneDistribute::New(mesh_);

    exchangeFields.setUpCommforZone(nextToInterface);

    Map<vector> mapCC =
        exchangeFields.getDatafromOtherProc(nextToInterface,mesh_.C());

    DynamicField<vector > cellCentre(100); // should be big enough avoids resizing

    const labelListList& CPCstencil = exchangeFields.getStencil();

    forAll(nextToInterface,celli)
    {
        if (nextToInterface[celli] && isCuboid_[celli])
        {
            cellCentre.clear();
            forAll(CPCstencil[celli],i)
            {
                const label& gblIdx = CPCstencil[celli][i];
                cellCentre.append
                (
                    exchangeFields.getValue(mesh_.C(),mapCC,gblIdx)
                );
            }

            if (CPCstencil[celli].size() == 27 || CPCstencil[celli].size() == 9)
            {
                sortStencilToIJKFormat
                (
                    celli,
                    CPCstencil[celli],
                    cellCentre,
                    mesh_.C()[celli]
                );

                // check if stencil was set correctly
                if (min(stencilHF_[celli]) == -1)
                {
                    stencilHF_[celli].setSize(0);
                    isCuboid_[celli] = false;
                }
            }
            else
            {
                stencilHF_[celli].setSize(0);
                isCuboid_[celli] = false;
            }

        }
        else
        {
            stencilHF_[celli].setSize(0);
        }

    }

}


Foam::Map< Foam::scalar > Foam::HFStencil::getDatafromOtherProc
(
    const boolList& nextToInterface,
    const volScalarField& alpha
)
{
    zoneDistribute& exchangeFields = zoneDistribute::New(mesh_);
    Map<scalar> MapAlpha=
        exchangeFields.getDatafromOtherProc(nextToInterface,alpha);
    return MapAlpha;
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


// ************************************************************************* //
