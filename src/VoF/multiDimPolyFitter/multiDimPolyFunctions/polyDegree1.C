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

#include "polyDegree1.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(polyDegree1, 0);
addToRunTimeSelectionTable(multiDimPolyFunctions, polyDegree1, word);

}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::polyDegree1::polyDegree1
(
    const Vector<label> dirs
)
: multiDimPolyFunctions(dirs)
{
    nTerms_ = 1;
    forAll(geomDir_,i)
    {
        if(geomDir_[i] == 1)
        {
            nTerms_ += 1;
        }
    }
    coeffs_.resize(nTerms_,0);
    termValues_.resize(nTerms_,0);
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


Foam::scalar Foam::polyDegree1::value(const vector& vec)
{
    // size has to be 4
    scalar value = coeffs_[0];
    forAll(geomDir_,i)
    {
        if(geomDir_[i] == 1)
        {
            value += coeffs_[i+1]*vec[i];
        }
    }

    return value;
}

const Foam::scalarField& Foam::polyDegree1::termValues(const vector& vec)
{
    termValues_[0]=1;
    label dimCounter = 0;
    forAll(geomDir_,i)
    {
        if(geomDir_[i] == 1)
        {
            dimCounter++;
            termValues_[dimCounter] = vec[i];
        }
    }

    return termValues_;
}

// ************************************************************************* //
