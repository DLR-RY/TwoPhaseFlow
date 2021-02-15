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

#include "zoneDistributePoints.H"
#include "DynamicField.H"
#include "syncTools.H"

// ************************************************************************* //

template <typename Type>
Type Foam::zoneDistributePoints::faceValue
(
    const GeometricField<Type, fvPatchField, volMesh>& phi,
    const label faceI
) const
{
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
Foam::Map<Foam::Field <Type> >
Foam::zoneDistributePoints::getPointDatafromOtherProc
(
    const boolList& zone,
    const GeometricField<Type, fvPatchField, volMesh>& phi
)
{
    const labelListList& cPoints = mesh_.pointCells();
    const labelListList& fPoints = mesh_.pointFaces();

    // cannot be on empty or wedge since they are ignored in syncGlobalIdxMap
    DynamicList<Type> valList(16);
    Map<Field <Type> > values;

    // Insert all cells adjacent to the point in stencil
    forAll(boundaryPoints_, i)
    {
        valList.clear();
        const label pI = boundaryPoints_[i];
        forAll(cPoints[pI], j)
        {
            const label& celli = cPoints[pI][j];
            valList.append(phi[celli]);
        }
        // insert non wedge empty proc patches in stencil
        forAll(fPoints[pI], fi)
        {
            const label& faceI = fPoints[pI][fi];
            if (!mesh_.isInternalFace(faceI))
            {
                if(validBoundaryFace_[faceI - mesh_.nInternalFaces()])
                {
                    valList.append(faceValue(phi,faceI));
                }
            }
        }
        const Field<Type> fieldVal(valList);

        values.insert(pI, fieldVal);
    }

    // syncNeiCellMap processor values are appended to the list
    syncTools::syncPointMap(mesh_, values, appendOp<Type>());
    return values;
}


// ************************************************************************* //
