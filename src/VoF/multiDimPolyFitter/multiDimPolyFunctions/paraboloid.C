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

#include "paraboloid.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(paraboloid, 0);
addToRunTimeSelectionTable(multiDimPolyFunctions, paraboloid, word);

}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::paraboloid::paraboloid
(
    const Vector<label> dirs
)
: multiDimPolyFunctions(dirs),
nDims_(0)
{

    nTerms_ = 0; // no constant values
    forAll(geomDir_,i)
    {
        if(geomDir_[i] == 1)
        {
            nTerms_ += 2; // without mixed terms e.g x+x^2
            nDims_++;
        }
    }
    // compute the mixed terms
    if(nDims_ == 1)
    {
        nTerms_ += 0; // 1 dim no mixed terms
    }
    else if(nDims_ == 2)
    {
        nTerms_ += 1; // 2 dim 1 mixed term
    }
    else
    {
        FatalErrorInFunction << "Only one and two dimensions are supported"
                    << abort(FatalError);
    }

    if((geomDir_[2] != -1) || geomDir_.size() != 3)
    {
        FatalErrorInFunction << "geomDir has to be (1 1 -1) or (1 -1 -1)"
                             << "geomDir is " << geomDir_
                    << abort(FatalError);

    }
    coeffs_.resize(nTerms_,0);
    termValues_.resize(nTerms_,0);

}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
// Foam::scalar Foam::paraboloid::value(const vector2D& vec)
// {
//     // size has to be 3
//     return coeffs_[0] + coeffs_[1]*vec.x() + coeffs_[2]*vec.y();
// }

Foam::scalar Foam::paraboloid::value(const vector& vec)
{
    // size has to be 4
    notImplemented("paraboloid::value is not implemented");
    return coeffs_[0] + coeffs_[1]*vec.x() + coeffs_[2]*vec.y() + coeffs_[3]*vec.y();
}

// Foam::scalarField Foam::paraboloid::termValues(const vector2D& vec)
// {
//     // size has to be 3
//     scalarField termVal(3,0);
//     termVal[0]=1;
//     termVal[1]=vec.x();
//     termVal[2]=vec.y();
//     return termVal;
// }

const Foam::scalarField& Foam::paraboloid::termValues(const vector& vec)
{
    // c0*x + c1*y + c2*x^2 + c3*y^2 + c4*x*y
    label dimCounter = 0;
    forAll(geomDir_,i)
    {
        if(geomDir_[i] == 1)
        {

            termValues_[dimCounter] = vec[i];
            termValues_[dimCounter + nDims_] = vec[i]*vec[i];
            dimCounter++;
        }
    }

    if(nDims_ == 2) // mixed terms
    {
        termValues_[4] = vec[0]*vec[1];
    }

    return termValues_;
}

// ************************************************************************* //
