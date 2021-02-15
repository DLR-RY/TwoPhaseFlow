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

#include "leastSquareInterpolate.H"

#include "emptyPolyPatch.H"
#include "processorPolyPatch.H"
#include "wedgePolyPatch.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
template<class T>
Foam::leastSquareInterpolate<T>::leastSquareInterpolate
(
    word functionName,
    Vector<label> geomDir
)
:
    polyFitter_(functionName,geomDir)
{

}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
template<class T>
T Foam::leastSquareInterpolate<T>::interpolate
(
    const List<vector>& positions,
    const List<T>& listValue
)
{
    List<T> fitData = polyFitter_.fitData
    (
        positions,
        listValue
    );

    return fitData[0];
}


template<class T>
Foam::Map < T >  Foam::leastSquareInterpolate<T>::interpolate
(
    const Map <List<vector> >& positions,
    const Map <List<T> >& listValue
)
{
    Map< T > interpolMap(positions.size());
    Map <List<vector> >::const_iterator iterPos = positions.cbegin();
    typename Map <List<T> >::const_iterator iterValue = listValue.cbegin();

    while(iterPos != positions.cend())
    {
        const List<vector>& positions = iterPos();
        const List<T>& listValue = iterValue();

        const T value = interpolate(positions,listValue);

        interpolMap.insert(iterPos.key(),value);

        ++iterPos;
        ++iterValue;

    }

    return interpolMap;
}

template class Foam::leastSquareInterpolate<Foam::scalar>;
template class Foam::leastSquareInterpolate<Foam::vector>;


// ************************************************************************* //
