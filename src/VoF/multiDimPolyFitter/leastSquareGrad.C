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

#include "leastSquareGrad.H"

#include "emptyPolyPatch.H"
#include "processorPolyPatch.H"
#include "wedgePolyPatch.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
template<class T>
Foam::leastSquareGrad<T>::leastSquareGrad
(
    word functionName,
    Vector<label> geomDir
)
:
    polyFitter_(functionName,geomDir),
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

}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
template<class T>
typename Foam::outerProduct<Foam::vector, T>::type Foam::leastSquareGrad<T>::grad
(
    const List<vector>& positions,
    const List<T>& listValue
)
{
    typedef typename outerProduct<vector, T>::type GradType;

    List<T> fitData = polyFitter_.fitData
    (
        positions,
        listValue
    );

    label dimCounter = 0;
    GradType grad = Zero;

    if(nDims_ == 3)
    {
        grad = GradType(fitData[1],fitData[2],fitData[3]);
    }
    else
    {
        forAll(geomDir_,i)
        {
            if(geomDir_[i] == 1)
            {
                grad[i] = fitData[dimCounter+1];
                dimCounter++;
            }
        }
    }
  

    return grad;
}

namespace Foam // needed g++ bug
{
    template<>
    tensor leastSquareGrad<vector>::grad
    (
        const List<vector>& positions,
        const List<vector>& listValue
    )
    {

        List<vector> fitData = polyFitter_.fitData
        (
            positions,
            listValue
        );

        label dimCounter = 0;
        tensor t = Zero;

        if(nDims_ == 3)
        {
            t = tensor(fitData[1],fitData[2],fitData[3]);
        }
        else
        {
            forAll(geomDir_,i)
            {
                if(geomDir_[i] == 1)
                {
                    const vector& fitVec = fitData[dimCounter+1];
                    forAll(fitVec,j)
                    {
                        t[i*3 + j] = fitVec[j];
                    }
                    dimCounter++;
                }
            }
        }


        return t;
    }
}

template<class T>
Foam::Map < typename Foam::outerProduct<Foam::vector, T>::type >  Foam::leastSquareGrad<T>::grad
(
    const Map <List<vector> >& positions,
    const Map <List<T> >& listValue
)
{
    typedef typename outerProduct<vector, T>::type GradType;

    Map< GradType > gradMap(positions.size());
    Map <List<vector> >::const_iterator iterPos = positions.cbegin();
    typename Map <List<T> >::const_iterator iterValue = listValue.cbegin();

    while(iterPos != positions.cend())
    {
        const List<vector>& positions = iterPos();
        const List<T>& listValue = iterValue();

        GradType grad = this->grad(positions,listValue);

        gradMap.insert(iterPos.key(),grad);

        ++iterPos;
        ++iterValue;

    }

    return gradMap;
}

template class Foam::leastSquareGrad<Foam::scalar>;
template class Foam::leastSquareGrad<Foam::vector>;


// ************************************************************************* //
