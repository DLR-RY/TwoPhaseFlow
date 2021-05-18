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

#include "surfaceTensionForceModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(surfaceTensionForceModel, 0);
    defineRunTimeSelectionTable(surfaceTensionForceModel, components);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::surfaceTensionForceModel::surfaceTensionForceModel
(
    const word& type,
    const dictionary& dict,
    const volScalarField& alpha1,
    const surfaceScalarField& phi,
    const volVectorField& U
)
:
    alpha1_(alpha1),
    phi_(phi),
    U_(U),
    surfTenModel_(nullptr),
    deltaFunctionModel_(nullptr),
    nHatf_
    (
        IOobject
        (
            "nHatf_",
            alpha1_.time().timeName(),
            alpha1_.mesh()
        ),
        alpha1_.mesh(),
        dimensionedScalar("nHatf", dimArea, 0.0)
    ),
    K_
    (
        IOobject
        (
            "K_",
            alpha1_.time().timeName(),
            alpha1_.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        alpha1_.mesh(),
        dimensionedScalar("K", dimless/dimLength, 0.0),
        "zeroGradient"
    ),
    Kf_
    (
        IOobject
        (
            "Kf_",
            alpha1_.time().timeName(),
            alpha1_.mesh()
        ),
        alpha1_.mesh(),
        dimensionedScalar("Kf", dimless/dimLength, 0.0)
    )

{

    surfTenModel_ = surfaceTensionModel::New(dict,alpha1.mesh());
    word deltaFunctionType =
        dict.lookupOrDefault<word>("deltaFunctionModel","alphaCSF");
    deltaFunctionModel_ = deltaFunctionModel::New(deltaFunctionType,dict,alpha1);
}





// ************************************************************************* //
