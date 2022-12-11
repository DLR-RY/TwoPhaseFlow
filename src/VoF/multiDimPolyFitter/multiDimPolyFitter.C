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

#include "multiDimPolyFitter.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class T>
Foam::multiDimPolyFitter<T>::multiDimPolyFitter
(
    const word polyFunctionName,
    const Vector<label> geomDirs
)
:
    polyFunc_(multiDimPolyFunctions::New(polyFunctionName,geomDirs)),
    A_(polyFunc_->nTerms(), 0.0, Zero)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class T>
void Foam::multiDimPolyFitter<T>::resetMatrix()
{
    for (label i=0;i<A_.size();i++)
    {
        A_.data()[i] = Zero;
    }
    A_.source() = Zero;
}


template<class T>
void Foam::multiDimPolyFitter<T>::fillMatrix
(
    const scalarField& polyTerms,
    const T& value
)
{
    label size = A_.n();

    //simple matrix is a square
    scalar tmpValue = 0;
    for(label i=0; i<size; ++i) // col
    {
        A_.source()[i] += polyTerms[i]*value;
        scalar* luMatrixi = A_[i];
        tmpValue = polyTerms[i];

        for(label j=0; j<size; ++j) // row
        {
            luMatrixi[j] += tmpValue*polyTerms[j];
        }
    }
}


template<class T>
void Foam::multiDimPolyFitter<T>::fillMatrix
(
    const scalarField& polyTerms,
    const T& value,
    const scalar weight
)
{
    label size = A_.n();

    //simple matrix is a square
    scalar tmpValue = 0;
    for (label i=0; i<size; ++i) // col
    {
        A_.source()[i] += polyTerms[i]*value*weight;
        scalar* __restrict luMatrixi = A_[i];
        tmpValue = polyTerms[i];

        for (label j=0; j<size; j++) // row
        {
            luMatrixi[j] += tmpValue*polyTerms[j]*weight;
        }
    }
}


template<class T>
void Foam::multiDimPolyFitter<T>::fillMatrix
(
    const scalarField& polyTerms,
    scalarSymmetricSquareMatrix& A
)
{
    label size = A.n();

    // simple matrix is a square
    for(label i=0; i<size; ++i)
    {
        for(label j=0; j<size; ++j)
        {
            A[i][j] += polyTerms[i]*polyTerms[j];
        }
    }
}


template<class T>
Foam::Field<T> Foam::multiDimPolyFitter<T>::fitData
(
    const List<scalarField>& listPolyTerms,
    const List<T>& listValue
)
{
    // operator= does not work
    resetMatrix();
    if (listPolyTerms.size() == listValue.size())
    {
        forAll(listPolyTerms,i)
        {
            fillMatrix
            (
                listPolyTerms[i],
                listValue[i]
            );
        }
        // Solve the matrix using Gaussian elimination with pivoting
        return A_.LUsolve();
    }
    else
    {
        FatalErrorInFunction
            << "size of listPolyTerms: " << listPolyTerms.size()
            << "size of listValues is:" <<  listValue.size()
            << "they have to match"
            << exit(FatalError);
        return Field<T>(0);
    }
}


template<class T>
Foam::Field<T> Foam::multiDimPolyFitter<T>::fitData
(
    const List<scalarField>& listPolyTerms,
    const List<T>& listValue,
    const List<scalar>& listWeight
)
{
    // operator= does not work
    resetMatrix();
    if (listPolyTerms.size() == listValue.size())
    {
        forAll(listPolyTerms, i)
        {
            fillMatrix
            (
                listPolyTerms[i],
                listValue[i],
                listWeight[i]
            );
        }

        // Solve the matrix using Gaussian elimination with pivoting
        return A_.LUsolve();
    }
    else
    {
        FatalErrorInFunction
            << "size of listPolyTerms: " << listPolyTerms.size()
            << "size of listValues is:" <<  listValue.size()
            << "they have to match"
            << exit(FatalError);
        return Field<T>(0);
    }
}


template<class T>
Foam::scalarSymmetricSquareMatrix Foam::multiDimPolyFitter<T>::computeInverse
(
    const List<scalarField>& listPolyTerms
)
{
    // operator= does not work
    scalarSymmetricSquareMatrix symMatrix(A_.n(), 0);
    forAll(listPolyTerms,i)
    {
        fillMatrix
        (
            listPolyTerms[i],
            symMatrix
        );
    }

    return inv(symMatrix);
}


template<class T>
Foam::Field<T> Foam::multiDimPolyFitter<T>::computeMatrixSource
(
    const List<scalarField>& listPolyTerms,
    const List<T>& listValue
)
{
    if (listPolyTerms.size() != listValue.size())
    {
        FatalErrorInFunction
            << "size of listPolyTerms: " << listPolyTerms.size()
            << "size of listValues is:" <<  listValue.size()
            << "they have to match"
            << exit(FatalError);
    }

    Field<T> source(listPolyTerms.size(), Zero);

    forAll(source, i)
    {
        forAll(listPolyTerms[i], j)
        {
            source[i] +=listPolyTerms[i][j]*listValue[i];
        }
    }

    return source;
}


template<class T>
Foam::Field<T> Foam::multiDimPolyFitter<T>::fitData
(
    const List<vector>& positions,
    const List<T>& listValue
)
{
    // operator= does not work
    if (positions.size() != listValue.size())
    {
        FatalErrorInFunction
            << "size of positions and listValues don't match" << nl
            << "size of positions is: " << positions.size()  << nl
            << "size of listValues is: " <<  listValue.size() << nl
            << exit(FatalError);
    }

    resetMatrix();

    forAll(positions, i)
    {
        fillMatrix
        (
            polyFunc_->termValues(positions[i]),
            listValue[i]
        );
    }

    // Solve the matrix using Gaussian elimination with pivoting
    return A_.LUsolve();

}


template<class T>
Foam::Field<T> Foam::multiDimPolyFitter<T>::fitData
(
    const List<vector>& positions,
    const List<T>& listValue,
    const List<scalar>& listWeight
)
{
    // operator= does not work
    if (positions.size() != listValue.size())
    {
        FatalErrorInFunction
            << "size of positions and listValues don't match" << nl
            << "size of positions is: " << positions.size()  << nl
            << "size of listValues is: " <<  listValue.size() << nl
            << exit(FatalError);
    }

    resetMatrix();

    forAll(positions, i)
    {
        fillMatrix
        (
            polyFunc_->termValues(positions[i]),
            listValue[i],
            listWeight[i]
        );
    }

    // Solve the matrix using Gaussian elimination with pivoting
    return A_.LUsolve();
}


template<class T>
Foam::scalarSymmetricSquareMatrix Foam::multiDimPolyFitter<T>::computeInverse
(
    const List<vector>& positions
)
{
    // operator= does not work
    scalarSymmetricSquareMatrix symMatrix(A_.n(), 0);
    forAll(positions, i)
    {
        fillMatrix
        (
            polyFunc_->termValues(positions[i]),
            symMatrix
        );
    }

    return inv(symMatrix);
}


template<class T>
Foam::Field<T> Foam::multiDimPolyFitter<T>::computeMatrixSource
(
    const List<vector>& positions,
    const List<T>& listValue
)
{
    if (positions.size() != listValue.size())
    {
        FatalErrorInFunction
            << "size of positions: " << positions.size()
            << "size of listValues is:" <<  listValue.size()
            << "they have to match"
            << exit(FatalError);
    }

    Field<T> source(polyFunc_->nTerms(), Zero);

    forAll(source, i)
    {
        scalarField polyTerms = polyFunc_->termValues(positions[i]);
        forAll(polyTerms, j)
        {
            source[i] += polyTerms[j]*listValue[i];
        }
    }

    return source;
}


template class Foam::multiDimPolyFitter<Foam::scalar>;
template class Foam::multiDimPolyFitter<Foam::vector>;

// ************************************************************************* //
