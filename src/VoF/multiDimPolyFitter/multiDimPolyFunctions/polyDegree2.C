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

#include "polyDegree2.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(polyDegree2, 0);
addToRunTimeSelectionTable(multiDimPolyFunctions, polyDegree2, word);

}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::polyDegree2::polyDegree2
(
    const Vector<label> dirs
)
: multiDimPolyFunctions(dirs),
nDims_(0)
{

    nTerms_ = 1;
    forAll(geomDir_,i)
    {
        if(geomDir_[i] == 1)
        {
            nTerms_ += 2; // without mixed terms
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
    else if(nDims_ == 3)
    {
        nTerms_ += 3; // 3 dim 3 mixed terms
    }
    coeffs_.resize(nTerms_,0);
    termValues_.resize(nTerms_,0);

}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
// Foam::scalar Foam::polyDegree2::value(const vector2D& vec)
// {
//     // size has to be 3
//     return coeffs_[0] + coeffs_[1]*vec.x() + coeffs_[2]*vec.y();
// }

Foam::scalar Foam::polyDegree2::value(const vector& vec)
{
    // size has to be 4
    notImplemented("polyDegree2::value is not implemented");
    return coeffs_[0] + coeffs_[1]*vec.x() + coeffs_[2]*vec.y() + coeffs_[3]*vec.y();
}

// Foam::scalarField Foam::polyDegree2::termValues(const vector2D& vec)
// {
//     // size has to be 3
//     scalarField termVal(3,0);
//     termVal[0]=1;
//     termVal[1]=vec.x();
//     termVal[2]=vec.y();
//     return termVal;
// }

const Foam::scalarField& Foam::polyDegree2::termValues(const vector& vec)
{
    // c0 + c1*x + c2*y + c3*z + c4*x^2 + c5*y^2 + c6*z^2 + c7*x*y  + c8*x*z + c9*y*z
    termValues_[0]=1;
    label dimCounter = 0;
    forAll(geomDir_,i)
    {
        if(geomDir_[i] == 1)
        {
            dimCounter++;
            termValues_[dimCounter] = vec[i];
            termValues_[dimCounter + nDims_] = vec[i]*vec[i];
        }
    }

    // account for mixed terms
    // 1 dim do nothing
    if(nDims_ == 3)
    {
        termValues_[7] = vec.x()*vec.y();
        termValues_[8] = vec.x()*vec.z();
        termValues_[9] = vec.y()*vec.z();
    }
    else if(nDims_ == 2)
    {

        bool firstDir = true;
        forAll(geomDir_,i)
        {
            if(geomDir_[i] == 1)
            {
                if(firstDir)
                {
                    termValues_[5] = vec[i];
                    firstDir = false;
                }
                else
                {
                    termValues_[5] *= vec[i];
                }
            }
        }
    }

    return termValues_;
}

// ************************************************************************* //
