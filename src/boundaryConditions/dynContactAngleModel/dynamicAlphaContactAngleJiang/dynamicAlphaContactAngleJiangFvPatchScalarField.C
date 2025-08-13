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

#include "dynamicAlphaContactAngleJiangFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volMesh.H"

#include "volFields.H"
#include "immiscibleIncompressibleTwoPhaseMixture.H"
#include "incompressibleTwoPhaseMixture.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::dynamicAlphaContactAngleJiangFvPatchScalarField::
dynamicAlphaContactAngleJiangFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    alphaContactAngleTwoPhaseFvPatchScalarField(p, iF),
    theta0_(0.0)
{}


Foam::dynamicAlphaContactAngleJiangFvPatchScalarField::
dynamicAlphaContactAngleJiangFvPatchScalarField
(
    const dynamicAlphaContactAngleJiangFvPatchScalarField& gcpsf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    alphaContactAngleTwoPhaseFvPatchScalarField(gcpsf, p, iF, mapper),
    theta0_(gcpsf.theta0_)
{}


Foam::dynamicAlphaContactAngleJiangFvPatchScalarField::
dynamicAlphaContactAngleJiangFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    alphaContactAngleTwoPhaseFvPatchScalarField(p, iF, dict),
    theta0_(readScalar(dict.lookup("theta0")))
{
    evaluate();
}


Foam::dynamicAlphaContactAngleJiangFvPatchScalarField::
dynamicAlphaContactAngleJiangFvPatchScalarField
(
    const dynamicAlphaContactAngleJiangFvPatchScalarField& gcpsf
)
:
    alphaContactAngleTwoPhaseFvPatchScalarField(gcpsf),
    theta0_(gcpsf.theta0_)
{}


Foam::dynamicAlphaContactAngleJiangFvPatchScalarField::
dynamicAlphaContactAngleJiangFvPatchScalarField
(
    const dynamicAlphaContactAngleJiangFvPatchScalarField& gcpsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    alphaContactAngleTwoPhaseFvPatchScalarField(gcpsf, iF),
    theta0_(gcpsf.theta0_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::scalarField>
Foam::dynamicAlphaContactAngleJiangFvPatchScalarField::theta
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

    // muss selber erstellt werden weil rho noch nicht registriert ist
    const volScalarField& nu1 =
        this->db().objectRegistry::lookupObject<volScalarField>("nu1");

    const dictionary& transportProperties =
    this->db().objectRegistry::lookupObject<IOdictionary>
    (
        "transportProperties"
    );


    word phase1Name (wordList(transportProperties.lookup("phases"))[0]);
    word phase2Name (wordList(transportProperties.lookup("phases"))[1]);

    dimensionedScalar rho1 =
        transportProperties.subDict(phase1Name).get<dimensionedScalar>("rho");

    dimensionedScalar sigmap =
        transportProperties.get<dimensionedScalar>("sigma");

    const fvPatchScalarField& nu1p = nu1.boundaryField()[patchi];

    scalarField Ca(nu1p*rho1.value()*uwall/sigmap.value());

    if(this->db().objectRegistry::foundObject<volScalarField>("contactAngle"))
    {
        volScalarField &contactAngle=const_cast<volScalarField&>(
        this->db().objectRegistry::lookupObject<volScalarField>("contactAngle"));

        contactAngle.boundaryFieldRef()[patchi] =
        min
        (
            180/constant::mathematical::pi*acos(cos(theta0_)
            - tanh(scalar(4.96)*pow(pos(Ca)*Ca,0.706))
            *(cos(theta0_)+1)),
            scalar(180)
        );
    }

    if(this->db().objectRegistry::foundObject<volScalarField>("Ca"))
    {
        volScalarField &CaF=const_cast<volScalarField&>(
        this->db().objectRegistry::lookupObject<volScalarField>("Ca"));

        CaF.boundaryFieldRef()[patchi] = Ca;
    }

    return min
    (
        180/constant::mathematical::pi*acos(cos(theta0_)
        -tanh(scalar(4.96)*pow(pos(Ca)*Ca,0.706))
        *(cos(theta0_)+1)),
        scalar(180)
    );
}


void Foam::dynamicAlphaContactAngleJiangFvPatchScalarField::write(Ostream& os) const
{
    alphaContactAngleTwoPhaseFvPatchScalarField::write(os);
    os.writeKeyword("theta0") << theta0_ << token::END_STATEMENT << nl;
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        dynamicAlphaContactAngleJiangFvPatchScalarField
    );
}


// ************************************************************************* //
