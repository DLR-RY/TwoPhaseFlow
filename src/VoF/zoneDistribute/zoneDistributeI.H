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

#include "zoneDistribute.H"
#include "DynamicField.H"
#include "syncTools.H"

// ************************************************************************* //

template <typename Type>
Type Foam::zoneDistribute::getLocalValue
(
    const GeometricField<Type, fvPatchField, volMesh>& phi,
    const label& localIdx
) const
{
    if(localIdx < mesh_.nCells()) // internal: cellI
    {
        return phi[localIdx];
    }
    else
    {
        return faceValue(phi,localIdx);
    }

}

template <typename Type>
Type Foam::zoneDistribute::faceValue
(
    const GeometricField<Type, fvPatchField, volMesh>& phi,
    const label& localIdx
) const
{
    const label faceI =  localIdx + mesh_.nInternalFaces() - mesh_.nCells();

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

    return phi.boundaryField()[patchI][patchFaceI];
}

template <typename Type>
Type Foam::zoneDistribute::getValue
(
    const GeometricField<Type, fvPatchField, volMesh>& phi,
    const Map<Type>& valuesFromOtherProc,
    const label& gblIdx
) const
{
    if(gblIdx_.isLocal(gblIdx))
    {
        const label idx = gblIdx_.toLocal(gblIdx);
        return getLocalValue(phi,idx);
    }
    else // from other proc
    {
        return valuesFromOtherProc[gblIdx];
    }
}



template <typename Type>
Foam::Map<Foam::Field <Type> > Foam::zoneDistribute::getFields
(
    const boolList& zone,
    const GeometricField<Type, fvPatchField, volMesh>& phi
)
{
    zoneCPCStencil& stencil = zoneCPCStencil::New(mesh_);
    if(zone.size() != phi.size())
    {
        FatalErrorInFunction
        << "size of zone: " << zone.size()
        << "size of phi:" <<  phi.size()
        << "do not match. Did the mesh change?"
        << exit(FatalError);
        return Map<Field <Type> >(0);
    }

    // get values from other proc
    Map<Type> neiValues = getDatafromOtherProc
    (
        zone,
        phi
    );

    Map< Field <Type> > stencilWithValues;

    DynamicField< Type > tmpField(100);

    forAll(zone,celli)
    {
        if(zone[celli])
        {
            tmpField.clearStorage();

            forAll(stencil[celli],i)
            {
                const label& gblIdx = stencil[celli][i];
                tmpField.append(getValue(phi,neiValues,gblIdx));
            }
            stencilWithValues.insert(celli,tmpField);
        }

    }

    return stencilWithValues;

}

template <typename Type>
Foam::Map <Type>  Foam::zoneDistribute::getDatafromOtherProc
(
    const boolList& zone,
    const GeometricField<Type, fvPatchField, volMesh>& phi
)
{
    if(zone.size() != phi.size())
    {
        FatalErrorInFunction
        << "size of zone: " << zone.size()
        << "size of phi:" <<  phi.size()
        << "do not match. Did the mesh change?"
        << exit(FatalError);
        return Map<Type> (0);
    }

    // get values from other proc
    Map<Type> neiValues;
    List < Map < Type > > sendValues(Pstream::nProcs());

    if (Pstream::parRun())
    {
        forAll(send_,domaini)
        {
            const labelHashSet& toNei = send_[domaini];

            forAllConstIters(toNei,iter)
            {
                sendValues[domaini].insert
                (
                    iter.key(),
                    getLocalValue(phi,gblIdx_.toLocal(iter.key()))
                );

            }
        }

        PstreamBuffers pBufs(Pstream::commsTypes::nonBlocking);

        // Stream data into buffer
        for (label domain = 0; domain < Pstream::nProcs(); domain++)
        {

            if (domain != Pstream::myProcNo())
            {
                // Put data into send buffer
                UOPstream toDomain(domain, pBufs);

                toDomain << sendValues[domain];
            }
        }

        // wait until everything is written.
        pBufs.finishedSends();

        Map<Type> tmpValue;

        for (label domain = 0; domain < Pstream::nProcs(); domain++)
        {
            tmpValue.clear();

            if (domain != Pstream::myProcNo())
            {
                // get data from send buffer
                UIPstream fromDomain(domain, pBufs);

                fromDomain >> tmpValue;

                forAllConstIters(tmpValue, iter)
                {
                    neiValues.insert(iter.key(), iter());
                }
            }
        }
    }

    return neiValues;

}



// ************************************************************************* //
