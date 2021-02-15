/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2019 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
                            | Copyright (C) 2019 DLR
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

Class
    Foam::multiDimPolyFitter

Description
    Fit a polynominal function with the help of multiDimPolyFunction

    Original code supplied by Henning Scheufler, DLR (2019)

SourceFiles
    multiDimPolyFitter.C

\*---------------------------------------------------------------------------*/

#ifndef multiDimPolyFitter_H
#define multiDimPolyFitter_H

#include "simpleMatrix.H"
#include "multiDimPolyFunctions.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class multiDimPolyFitter Declaration
\*---------------------------------------------------------------------------*/

template<class T>
class multiDimPolyFitter
{
    // Private data

        //- Mesh whose cells and faces to cut at their intersection with an
        //-  isosurface.
        autoPtr<multiDimPolyFunctions> polyFunc_;

        simpleMatrix<T> A_;


public:

    //- Construct from number of terms
    multiDimPolyFitter
    (
        const word polyFunctionName,
        const Vector<label> geomDirs
    );

    void resetMatrix();

    label nCoeffs() const
    {
        return A_.n();
    }

    void fillMatrix
    (
        const scalarField& polyTerms,
        const T& value
    );

    void fillMatrix
    (
        const scalarField& polyTerms,
        const T& value,
        const scalar weight
    );

    void fillMatrix
    (
        const scalarField& polyTerms,
        scalarSymmetricSquareMatrix& A
    );

    //- Fit data
    Field<T> fitData
    (
        const List<scalarField>& listPolyTerms,
        const List<T>& listValue
    );

    //- Fit data
    Field<T> fitData
    (
        const List<scalarField>& listPolyTerms,
        const List<T>& listValue,
        const List<scalar>& listWeight
    );

    //- Compute inverse
    scalarSymmetricSquareMatrix computeInverse
    (
        const List<scalarField>& listPolyTerms
    );

    //- Compute source
    Field<T> computeMatrixSource
    (
        const List<scalarField>& listPolyTerms,
        const List<T>& listValue
    );

    //- Fit data
    Field<T> fitData
    (
        const List<vector>& positions,
        const List<T>& listValue
    );

    //- Fit data
    Field<T> fitData
    (
        const List<vector>& positions,
        const List<T>& listValue,
        const List<scalar>& listWeight
    );

    //- Compute inverse
    scalarSymmetricSquareMatrix computeInverse
    (
        const List<vector>& positions
    );

    //- Compute source
    Field<T> computeMatrixSource
    (
        const List<vector>& positions,
        const List<T>& listValue
    );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
