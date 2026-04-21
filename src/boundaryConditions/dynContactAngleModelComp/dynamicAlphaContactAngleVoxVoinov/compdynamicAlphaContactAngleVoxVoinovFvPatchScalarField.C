/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
     \\/     M anipulation  |
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

#include "compdynamicAlphaContactAngleVoxVoinovFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volMesh.H"


#include "volFields.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::compdynamicAlphaContactAngleVoxVoinovFvPatchScalarField::
compdynamicAlphaContactAngleVoxVoinovFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    alphaContactAngleTwoPhaseFvPatchScalarField(p, iF),
    init_(true),
    theta0_(0.0),
    ct_(0.0),
    muName_("")
{}


Foam::compdynamicAlphaContactAngleVoxVoinovFvPatchScalarField::
compdynamicAlphaContactAngleVoxVoinovFvPatchScalarField
(
    const compdynamicAlphaContactAngleVoxVoinovFvPatchScalarField& gcpsf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    alphaContactAngleTwoPhaseFvPatchScalarField(gcpsf, p, iF, mapper),
    init_(true),
    theta0_(gcpsf.theta0_),
    ct_(gcpsf.ct_),
    muName_(gcpsf.muName_)
{}


Foam::compdynamicAlphaContactAngleVoxVoinovFvPatchScalarField::
compdynamicAlphaContactAngleVoxVoinovFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    alphaContactAngleTwoPhaseFvPatchScalarField(p, iF, dict),
    init_(true),
    theta0_(readScalar(dict.lookup("theta0"))),
    ct_(readScalar(dict.lookup("ct"))),
    muName_(dict.lookup("muFieldName"))
{
    evaluate();
}


Foam::compdynamicAlphaContactAngleVoxVoinovFvPatchScalarField::
compdynamicAlphaContactAngleVoxVoinovFvPatchScalarField
(
    const compdynamicAlphaContactAngleVoxVoinovFvPatchScalarField& gcpsf
)
:
    alphaContactAngleTwoPhaseFvPatchScalarField(gcpsf),
    init_(true),
    theta0_(gcpsf.theta0_),
    ct_(gcpsf.ct_),
    muName_(gcpsf.muName_)
{}


Foam::compdynamicAlphaContactAngleVoxVoinovFvPatchScalarField::
compdynamicAlphaContactAngleVoxVoinovFvPatchScalarField
(
    const compdynamicAlphaContactAngleVoxVoinovFvPatchScalarField& gcpsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    alphaContactAngleTwoPhaseFvPatchScalarField(gcpsf, iF),
    init_(true),
    theta0_(gcpsf.theta0_),
    ct_(gcpsf.ct_),
    muName_(gcpsf.muName_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::scalarField>
Foam::compdynamicAlphaContactAngleVoxVoinovFvPatchScalarField::theta
(
    const fvPatchVectorField& Up,
    const fvsPatchVectorField& nHat
) const
{
    const vectorField nf(patch().nf());

    // Calculated the component of the velocity parallel to the wall
    vectorField Uwall(Up.patchInternalField() - Up);
    Uwall -= (nf & Uwall)*nf;

    // Find the direction of the interface parallel to the wall
    vectorField nWall(nHat - (nf & nHat)*nf);

    // Normalise nWall
    nWall /= (mag(nWall) + SMALL);

    // Calculate Uwall resolved normal to the interface parallel to
    // the interface
    scalarField uwall(-nWall & Uwall);

    const label patchi = this->patch().index();

    const volScalarField& mu1 = this->db().objectRegistry::lookupObject<volScalarField>(muName_);

    const dictionary& thermophysicalProperties = this->db().objectRegistry::lookupObject<IOdictionary>
    (
        "thermophysicalProperties"
    );

    dimensionedScalar sigmap(thermophysicalProperties.get<dimensionedScalar>("sigma"));

    const fvPatchScalarField&  mu1p = mu1.boundaryField()[patchi];

    scalarField Ca(mu1p*uwall/sigmap.value());

    // thetaD = (ct*  Ca + theta0^3)^(1/3)

    //  Ca^1/3 in rad
    if(this->db().objectRegistry::foundObject<volScalarField>("contactAngle"))
    {
        volScalarField &contactAngle=const_cast<volScalarField&>(
        this->db().objectRegistry::lookupObject<volScalarField>("contactAngle"));

        contactAngle.boundaryFieldRef()[patchi] = min
        (
            180/constant::mathematical::pi
            *(
                pow(ct_*pos(Ca)*Ca
              + pow(theta0_*180/constant::mathematical::pi,3),0.3333333)
            ),
            scalar(180)
        );
    }

    if(this->db().objectRegistry::foundObject<volScalarField>("Ca"))
    {
        volScalarField &CaF=const_cast<volScalarField&>
        (
            this->db().objectRegistry::lookupObject<volScalarField>("Ca")
        );
        CaF.boundaryFieldRef()[patchi] = Ca;
    }

    return min
    (
        180/constant::mathematical::pi
        *(
            pow(ct_*pos(Ca)*Ca
          + pow(theta0_*180/constant::mathematical::pi,3),0.3333333)
        ),
        scalar(180)
    );
}


void Foam::compdynamicAlphaContactAngleVoxVoinovFvPatchScalarField::write(Ostream& os) const
{
    alphaContactAngleTwoPhaseFvPatchScalarField::write(os);
    os.writeKeyword("theta0") << theta0_ << token::END_STATEMENT << nl;
    os.writeKeyword("ct") << ct_ << token::END_STATEMENT << nl;
    os.writeKeyword("muFieldName") << muName_ << token::END_STATEMENT << nl;
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        compdynamicAlphaContactAngleVoxVoinovFvPatchScalarField
    );
}


// ************************************************************************* //
