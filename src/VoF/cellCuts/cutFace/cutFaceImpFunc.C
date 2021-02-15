/*---------------------------------------------------------------------------*\
    Modified work | Copyright (c) 2017-2019, German Aerospace Center (DLR)
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

#include "cutFaceImpFunc.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::scalar Foam::cutFaceImpFunc::bisection(point p0, point p1)
{
    // return
    vector start = p0;
    vector end = p1;
    scalar a = func_.value(p0);
    scalar b = func_.value(p1);
    if (a * b >= 0)
    {
        return GREAT;
    }

    scalar c = a;
    vector p3 = p0;
    while (mag(p1-p0) >= 1e-12) // SMALL 1e-15
    {
        // Find middle point
        p3 = 0.5*(p0+p1);
        c = func_.value(p3);

        // Check if middle point is root
        if (mag(c) <= 1e-12)
        {
            break;


        // Decide the side to repeat the steps
        }
        else if (c*a < 0)
        {
            b = c;
            p1 = p3;
        }
        else
        {
            a = c;
            p0 = p3;
        }
    }

    return mag(start-p3)/mag(start-end);


}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::cutFaceImpFunc::cutFaceImpFunc(const fvMesh& mesh, scalarField& f, implicitFunction& func)
:
    cutFace(mesh),
    mesh_(mesh),
    f_(f),
    func_(func),
    subFaceCentre_(vector::zero),
    subFaceArea_(vector::zero),
    subFacePoints_(10),
    surfacePoints_(4),
    pointStatus_(10),
    weight_(10),
    faceStatus_(-1)
{
    clearStorage();
}

// * * * * * * * * * * * Public Member Functions  * * * * * * * * * * * * * //

Foam::label Foam::cutFaceImpFunc::calcSubFace
(
    const label& faceI,
    const scalar& cutValue
)
{
    clearStorage();
    const face& f = mesh_.faces()[faceI];
    label inLiquid = 0;
    label firstFullySubmergedPoint = -1;

    // loop over face
    forAll(f, i)
    {
        // pointStatus is cutValue - f_[f[i]] is the other way around in cutFaceIso
        pointStatus_.append( f_[f[i]] - cutValue );

        if (mag(pointStatus_[i]) < 10 * SMALL)
        {
            pointStatus_[i] = 0;
        }

        if (pointStatus_[i] > 10 * SMALL)
        {
            inLiquid++;
            if (firstFullySubmergedPoint == -1)
            {
                firstFullySubmergedPoint = i;
            }
        }
    }

    // calculate weight
    forAll(f, i)
    {

        const label nextIdx = (i+1) % f.size();
        const label pI = f[i];
        const label nextpI = f[nextIdx];
        if
        (
            (pointStatus_[i] < 0 && pointStatus_[nextIdx] > 0) ||
            (pointStatus_[i] > 0 && pointStatus_[nextIdx] < 0)
        ) // cut on edge cut Value is zero
        {
            scalar val = bisection(mesh_.points()[pI],mesh_.points()[nextpI]);
            weight_.append(val);
        }
        else
        {
            weight_.append(GREAT);
        }
    }


    if (inLiquid == f.size()) // fluid face
    {
        faceStatus_ = -1;
        subFaceCentre_ = mesh_.faceCentres()[faceI];
        subFaceArea_ = mesh_.faceAreas()[faceI];
        return faceStatus_;
    }
    else if (inLiquid == 0) // gas face
    {
        faceStatus_ = 1;
        subFaceCentre_ = vector::zero;
        subFaceArea_ = vector::zero;
        return faceStatus_;
    }


    cutFace::calcSubFace
    (
        faceI,
        pointStatus_,
        weight_,
        firstFullySubmergedPoint,
        subFacePoints_,
        surfacePoints_,
        faceStatus_,
        subFaceCentre_,
        subFaceArea_
    );

    return faceStatus_;
}

Foam::point Foam::cutFaceImpFunc::subFaceCentre()
{
    return subFaceCentre_;
}

Foam::vector Foam::cutFaceImpFunc::subFaceArea()
{
    return subFaceArea_;
}

Foam::DynamicList<Foam::point>& Foam::cutFaceImpFunc::subFacePoints()
{
    return subFacePoints_;
}

Foam::DynamicList<Foam::point>& Foam::cutFaceImpFunc::surfacePoints()
{
    return surfacePoints_;
}

void Foam::cutFaceImpFunc::clearStorage()
{
    subFaceCentre_ = vector::zero;
    subFaceArea_ = vector::zero;
    subFacePoints_.clear();
    surfacePoints_.clear();
    pointStatus_.clear();
    weight_.clear();
    faceStatus_ = -1;
}

// ************************************************************************* //
