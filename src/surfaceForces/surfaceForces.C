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

#include "surfaceForces.H"

#include "surfaceInterpolate.H"
#include "fvc.H"


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::surfaceForces::surfaceForces
(
    const volScalarField& alpha1,
    const surfaceScalarField& phi,
    const volVectorField& U,
    const dictionary& dict
)
:
    surfaceForcesCoeffs_(dict.subDict("surfaceForces")),
    alpha1_(alpha1),
    mesh_(alpha1.mesh()),
    surfTenForceModel_(nullptr),
    accModel_(nullptr)
{
    surfTenForceModel_ = surfaceTensionForceModel::New(surfaceForcesCoeffs_,alpha1,phi,U);
    accModel_ = accelerationForceModel::New(surfaceForcesCoeffs_,alpha1.mesh());
}



// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::surfaceScalarField> Foam::surfaceForces::surfaceTensionForce()
{
    return surfTenForceModel_->surfaceTensionForce();
}

Foam::tmp<Foam::surfaceScalarField> Foam::surfaceForces::accelerationForce()
{
    return accModel_->accelerationForce();
}

Foam::scalar Foam::surfaceForces::capillaryDt
(
    const volScalarField& rho1,
    const volScalarField& rho2
)
{
    volScalarField interface(pos0(alpha1_ - 0.001)*pos0(0.999 - alpha1_));
    volScalarField deltaX
    (
        interface*fvc::average(mag(mesh_.delta()))
      + neg0(interface-0.5)*dimensionedScalar("0",dimLength,GREAT)
    );

    dimensionedScalar smallSigma("smallSigma",dimensionSet(1,0,-2, 0, 0),SMALL);

    dimensionedScalar dt =
        gMin(sqrt((rho1 +rho2)*pow(deltaX.internalField(),3) /
        (2*constant::mathematical::pi*(sigma()().internalField()+smallSigma))));

    return dt.value();
}

Foam::scalar Foam::surfaceForces::capillaryDt
(
    const dimensionedScalar rho1,
    const dimensionedScalar rho2
)
{
    volScalarField  interface(pos0(alpha1_ - 0.001)*pos0(0.999 - alpha1_));
    volScalarField deltaX
    (
        interface*fvc::average(mag(mesh_.delta()))
      + neg0(interface-0.5)*dimensionedScalar("0",dimLength,GREAT)
    );
    dimensionedScalar smallSigma("smallSigma",dimensionSet(1,0,-2,0,0),SMALL);

    dimensionedScalar dt =
        min(sqrt((rho1 +rho2)*pow(deltaX,3) /
        (2*constant::mathematical::pi*(sigma()+smallSigma))));

    reduce(dt.value(),minOp<scalar>());
    return dt.value();
}

// ************************************************************************* //
