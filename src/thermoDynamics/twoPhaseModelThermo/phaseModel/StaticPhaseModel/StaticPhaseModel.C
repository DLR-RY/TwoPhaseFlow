/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2017 OpenCFD Ltd.
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

#include "StaticPhaseModel.H"


#include "fvcDdt.H"
#include "fvcDiv.H"
#include "surfaceInterpolate.H"


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasePhaseModel>
Foam::StaticPhaseModel<BasePhaseModel>::StaticPhaseModel
(
    const fvMesh& mesh,
    const dictionary& dict,
    const word& phaseName
)
:
    BasePhaseModel(mesh,dict, phaseName),
    U_(mesh.lookupObject<volVectorField>("U")),
    phi_
    (
        IOobject
        (
            IOobject::groupName("phi", phaseModel::name()),
            mesh.time().timeName(),
            mesh
        ),
        mesh,
        dimensionedScalar("0", dimensionSet(0, 3, -1, 0, 0), 0)
    ),
    alphaPhi_
    (
        IOobject
        (
            IOobject::groupName("alphaPhi", phaseModel::name()),
            mesh.time().timeName(),
            mesh
        ),
        mesh,
        dimensionedScalar("0", dimensionSet(0, 3, -1, 0, 0), 0)
    )
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasePhaseModel>
void Foam::StaticPhaseModel<BasePhaseModel>::correct()
{
    BasePhaseModel::correct();
}


template<class BasePhaseModel>
Foam::tmp<Foam::surfaceScalarField>
Foam::StaticPhaseModel<BasePhaseModel>::phi() const
{
    return tmp<surfaceScalarField>::New
    (
        IOobject
        (
            IOobject::groupName("phi", phaseModel::name()),
            U_.mesh().time().timeName(),
            U_.mesh()
        ),
        U_.mesh(),
        dimensionedScalar("0", dimensionSet(0, 3, -1, 0, 0), 0)
    );
}


template<class BasePhaseModel>
const Foam::surfaceScalarField&
Foam::StaticPhaseModel<BasePhaseModel>::phi()
{
    phi_ = dimensionedScalar("0", dimensionSet(0, 3, -1, 0, 0), 0);
    return phi_;
}


template<class BasePhaseModel>
Foam::tmp<Foam::surfaceScalarField>
Foam::StaticPhaseModel<BasePhaseModel>::alphaPhi() const
{
    return tmp<surfaceScalarField>
    (
        new surfaceScalarField
        (
            IOobject
            (
                IOobject::groupName("alphaPhi", phaseModel::name()),
                U_.mesh().time().timeName(),
                U_.mesh()
            ),
            U_.mesh(),
            dimensionedScalar("0", dimensionSet(0, 3, -1, 0, 0), 0)
        )
    );
}


template<class BasePhaseModel>
Foam::surfaceScalarField&
Foam::StaticPhaseModel<BasePhaseModel>::alphaPhi()
{
    alphaPhi_ = dimensionedScalar("0", dimensionSet(0, 3, -1, 0, 0), 0);
    return alphaPhi_;
}


template<class BasePhaseModel>
Foam::tmp<Foam::volVectorField>
Foam::StaticPhaseModel<BasePhaseModel>::U() const
{
    return tmp<volVectorField>
    (
        new volVectorField
        (
            IOobject
            (
                IOobject::groupName("U", phaseModel::name()),
                U_.mesh().time().timeName(),
                U_.mesh()
            ),
            U_.mesh(),
            dimensionedVector("zero", dimVelocity, vector::zero)
        )
    );
}


template<class BasePhaseModel>
Foam::tmp<Foam::surfaceScalarField> Foam::StaticPhaseModel<BasePhaseModel>
::diffNo() const
{
    tmp<surfaceScalarField> tkapparhoCpbyDelta
    (
        sqr(U_.mesh().surfaceInterpolation::deltaCoeffs())
       *fvc::interpolate(this->kappa().ref())
       /fvc::interpolate((this->Cp()*this->rho())())
    );

    return tkapparhoCpbyDelta;
}


// ************************************************************************* //
