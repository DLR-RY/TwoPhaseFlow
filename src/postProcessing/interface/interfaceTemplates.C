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

#include "interface.H"
#include "polyMesh.H"
#include "syncTools.H"
#include "surfaceFields.H"
#include "OFstream.H"
#include "meshTools.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //



template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::interface::interpolate
(
    const GeometricField<Type, fvPatchField, volMesh>& cCoords,
    const Field<Type>& pCoords
) const
{

    return Field<Type>(this->points().size()); //,Type::zero);
}

//template< >
//Foam::tmp<Foam::Field<Foam::scalar>>
//Foam::interface::interpolate
//(
//    const GeometricField<scalar, fvPatchField, volMesh>& cCoords,
//    const Field<scalar>& pCoords
//) const
//{

//    return Field<scalar>(this->points().size(),0);
//}

// ************************************************************************* //
