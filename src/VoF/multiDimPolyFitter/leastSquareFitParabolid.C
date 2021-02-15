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

#include "leastSquareFitParabolid.H"

#include "emptyPolyPatch.H"
#include "processorPolyPatch.H"
#include "wedgePolyPatch.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


Foam::cartesianCS Foam::leastSquareFitParabolid::makeLocalCoordSystem
(
    const point& basePoint,
    const vector& normal
)
{
    // find in plane vector

    vector e0(0,0,0);
    vector n = normal/mag(normal);

    if(nDims_ == 3)
    {
        vector base(1,0,0);

        scalar nComp = n & base;

        if (mag(nComp) > 0.8)
        {
            // Was bad guess. Try with different vector.
            base.x() = 0;
            base.y() = 1;

            nComp = n & base;

            if (mag(nComp) > 0.8)
            {
                base.y() = 0;
                base.z() = 1;

                nComp = n & base;
            }
        }

        // Use component normal to n as base vector.
        e0 = base - nComp*n;

        e0 /= mag(e0) + VSMALL;
    }
    else if(nDims_ == 2)
    {
        // very stupid solution
        if (geomDir_.x() == -1)
        {
            e0 = n ^ vector(1,0,0);
        }
        else if (geomDir_.y() == -1)
        {
            e0 = n ^ vector(0,1,0);
        }
        else if (geomDir_.z() == -1)
        {
            e0 = n ^ vector(0,0,1);
        }

        e0 /= mag(e0) + VSMALL;
    }
    else
    {
        // Fatal Error
        FatalErrorInFunction
        << "the geometric dimension is not 2D or 2D" << nl
        << " the dimension is: " << nDims_
        << exit(FatalError);
    }

    cartesianCS pCS("planeCellCoord", basePoint, normal,e0);

    return pCS;


}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::leastSquareFitParabolid::leastSquareFitParabolid
(
    const Vector<label> geomDir,
    const Vector<label> explicitDim
)
:
    polyFitter_("paraboloid",explicitDim),
    geomDir_(geomDir),
    nDims_(0)
{
    // compute number of dimensions
    forAll(geomDir_,i)
    {
        if(geomDir_[i] == 1)
        {
            nDims_++;
        }
    }

    if(nDims_ <= 1)
    {
        // Fatal Error
        FatalErrorInFunction
        << "the geometric dimension is not 2D or 2D" << nl
        << " the dimension is: " << nDims_ << nl
        << exit(FatalError);
    }


}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
Foam::scalarField Foam::leastSquareFitParabolid::fitParaboloid
(
    const point& centre,
    const vector& normal,
    const vectorField& positions
)
{
    cartesianCS localCoord = makeLocalCoordSystem
    (
        centre,
        normal
    );
    vectorField localPositions (localCoord.localPosition(positions));

    List<scalar> listValue(localPositions.size());
    forAll(listValue,i)
    {
        listValue[i] = localPositions[i].z();
    }

    scalarField fitData = polyFitter_.fitData
    (
        localPositions,
        listValue
    );

    return fitData;
}

Foam::scalarField Foam::leastSquareFitParabolid::fitParaboloid
(
    const point& centre,
    const vector& normal,
    const vectorField& positions,
    const scalarField& weight
)
{
    cartesianCS localCoord = makeLocalCoordSystem
    (
        centre,
        normal
    );
    vectorField localPositions (localCoord.localPosition(positions));

    List<scalar> listValue(localPositions.size());
    forAll(listValue,i)
    {
        listValue[i] = localPositions[i].z();
    }

    scalarField fitData = polyFitter_.fitData
    (
        localPositions,
        listValue,
        weight
    );

    return fitData;
}

// Foam::Map < Foam::vector >  Foam::leastSquareFitParabolid::grad
// (
//     const Map <List<vector> >& positions,
//     const Map <List<scalar> >& listValue
// )
// {
//     Map< vector > gradMap(positions.size());
//     Map <List<vector> >::const_iterator iterPos = positions.cbegin();
//     Map <List<scalar> >::const_iterator iterValue = listValue.cbegin();

//     while(iterPos != positions.cend())
//     {
//         const List<vector>& positions = iterPos();
//         const List<scalar>& listValue = iterValue();

//         Info << "positions  " << positions <<endl;
//         Info << "listValue  " << listValue <<endl;
//         vector grad = this->grad(positions,listValue);

//         gradMap.insert(iterPos.key(),grad);

//         ++iterPos;
//         ++iterValue;

//     }

//     return gradMap;
// }


// ************************************************************************* //
