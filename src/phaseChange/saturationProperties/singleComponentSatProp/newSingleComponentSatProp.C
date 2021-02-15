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

#include "singleComponentSatProp.H"
//#include "twoPhaseModelThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::singleComponentSatProp>
Foam::singleComponentSatProp::New
(
    const fvMesh& mesh,
    const dictionary& dict

)
{

    word singleComponentSatPropTypeName
    (
        dict.get<word>("singleComponentSatProp")
    );

    Info<< "Selecting saturation model "
        << singleComponentSatPropTypeName << endl;

    componentsConstructorTable::iterator cstrIter =
        componentsConstructorTablePtr_
            ->find(singleComponentSatPropTypeName);

    if (cstrIter == componentsConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "saturationModel singleComponentSatProp::New"
        )   << "Unknown singleComponentSatProp type "
            << singleComponentSatPropTypeName << endl << endl
            << "Valid  singleComponentSatProps are : " << endl
            << componentsConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<singleComponentSatProp>(cstrIter()(mesh,dict));
}


// ************************************************************************* //
