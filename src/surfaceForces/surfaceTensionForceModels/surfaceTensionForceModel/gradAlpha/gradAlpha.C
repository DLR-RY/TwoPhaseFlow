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

#include "gradAlpha.H"
#include "addToRunTimeSelectionTable.H"

#include "alphaContactAngleTwoPhaseFvPatchScalarField.H"
#include "mathematicalConstants.H"
#include "surfaceInterpolate.H"
#include "fvcDiv.H"
#include "fvcGrad.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(gradAlpha, 0);
    addToRunTimeSelectionTable(surfaceTensionForceModel,gradAlpha, components);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::gradAlpha::gradAlpha
(
    const dictionary& dict,
    const volScalarField& alpha1,
    const surfaceScalarField& phi,
    const volVectorField& U
)
:
    surfaceTensionForceModel
    (
        typeName,
        dict,
        alpha1,
        phi,
        U
    ),
    deltaN_
    (
        "deltaN",
        1e-8/pow(average(alpha1.mesh().V()), 1.0/3.0)
    )
{

}


// * * * * * * * * * * * * * * Public Access Member Functions  * * * * * * * * * * * * * * //

void Foam::gradAlpha::correctContactAngle
(
    surfaceVectorField::Boundary& nHatb,
    surfaceVectorField::Boundary& gradAlphaf
)
{
    scalar convertToRad = Foam::constant::mathematical::pi/180.0;

    const fvMesh& mesh = alpha1_.mesh();
    const volScalarField::Boundary& abf = alpha1_.boundaryField();

    const fvBoundaryMesh& boundary = mesh.boundary();

    forAll(boundary, patchi)
    {
        if (isA<alphaContactAngleTwoPhaseFvPatchScalarField>(abf[patchi]))
        {
            alphaContactAngleTwoPhaseFvPatchScalarField& acap =
                const_cast<alphaContactAngleTwoPhaseFvPatchScalarField&>
                (
                    refCast<const alphaContactAngleTwoPhaseFvPatchScalarField>
                    (
                        abf[patchi]
                    )
                );

            fvsPatchVectorField& nHatp = nHatb[patchi];
            const scalarField theta
            (
                convertToRad*acap.theta(U_.boundaryField()[patchi], nHatp)
            );

            const vectorField nf
            (
                boundary[patchi].nf()
            );

            // Reset nHatp to correspond to the contact angle
            const scalarField a12(nHatp & nf);
            const scalarField b1(cos(theta));

            scalarField b2(nHatp.size());
            forAll(b2, facei)
            {
                b2[facei] = cos(acos(a12[facei]) - theta[facei]);
            }

            const scalarField det(1.0 - a12*a12);

            scalarField a((b1 - a12*b2)/det);
            scalarField b((b2 - a12*b1)/det);

            nHatp = a*nf + b*nHatp;
            nHatp /= (mag(nHatp) + deltaN_.value());

            acap.gradient() = (nf & nHatp)*mag(gradAlphaf[patchi]);
            acap.evaluate();
        }
    }
}


void Foam::gradAlpha::correct()
{
    deltaFunctionModel_->correct();

    const fvMesh& mesh = alpha1_.mesh();
    const surfaceVectorField& Sf = mesh.Sf();

    // Cell gradient of alpha
    const volVectorField gradAlpha(fvc::grad(alpha1_, "nHat"));

    // Interpolated face-gradient of alpha
    surfaceVectorField gradAlphaf(fvc::interpolate(gradAlpha));

    // Face unit interface normal
    surfaceVectorField nHatfv(gradAlphaf/(mag(gradAlphaf) + deltaN_));

    correctContactAngle(nHatfv.boundaryFieldRef(), gradAlphaf.boundaryFieldRef());

    // Face unit interface normal flux
    nHatf_ = nHatfv & Sf;

    // Simple expression for curvature
    K_ = -fvc::div(nHatf_);

    Kf_ = fvc::interpolate(K_);
}



// ************************************************************************* //
