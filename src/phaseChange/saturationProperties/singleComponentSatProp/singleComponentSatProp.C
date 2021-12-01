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
#include "zeroGradientFvPatchFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(singleComponentSatProp, 0);
    defineRunTimeSelectionTable(singleComponentSatProp, components);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::singleComponentSatProp::singleComponentSatProp
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    IOobject
    (
        "satProperties::singleComponentSatProp",
        mesh.time().constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    ),
    singleComponentSatPropCoeffs_(dict),
    mesh_(mesh),
    TSat_
    (
        IOobject
        (
            "TSat",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimTemperature,
        calculatedFvPatchScalarField::typeName
    ),
    pSat_
    (
        IOobject
        (
            "pSat",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimPressure,
        calculatedFvPatchScalarField::typeName
    ),
    L_
    (
        IOobject
        (
            "L",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimEnergy/dimMass,
        calculatedFvPatchScalarField::typeName
    ),
    Tmin_(modelDict().get<scalar>("Tmin")),
    Tmax_(modelDict().get<scalar>("Tmax"))
{
}

// * * * * * * * * * * * * * * Public Access Member Functions  * * * * * * * * * * * * * * //
void Foam::singleComponentSatProp::correct
(
    const volScalarField& p,
    const volScalarField& T
)
{
    // correct TSat
    forAll (TSat_.boundaryField(), patchi)
    {
        forAll (TSat_.boundaryField()[patchi], i)
        {
            TSat_.boundaryFieldRef()[patchi][i] =
            satT
            (
                p.boundaryField()[patchi][i]
            );
        }
    }

    forAll (TSat_, celli)
    {
        TSat_[celli] = satT(p[celli]);
    }


    // correct pSat
    forAll (pSat_.boundaryField(), patchi)
    {
        forAll (pSat_.boundaryField()[patchi], i)
        {
            pSat_.boundaryFieldRef()[patchi][i] =
            satP
            (
                T.boundaryField()[patchi][i]
            );
        }
    }

    forAll (pSat_, celli)
    {
        pSat_[celli] = satP(T[celli]);
    }

    // correct L
    forAll (L_.boundaryField(), patchi)
    {
        forAll (L_.boundaryField()[patchi], i)
        {
            L_.boundaryFieldRef()[patchi][i] =
            satL
            (
                p.boundaryField()[patchi][i]
            );
        }
    }

    forAll (L_, celli)
    {
        L_[celli] = satL(p[celli]);
    }
}

Foam::scalar
Foam::singleComponentSatProp::satP(scalar t) const
{
    notImplemented("bool Foam::singleComponentSatProp::satP(scalar t)");
    return 0.0;
}

Foam::scalar
Foam::singleComponentSatProp::satT(scalar p) const
{
    notImplemented("bool Foam::singleComponentSatProp::satT(scalar p)");
    return 0.0;
}

Foam::scalar
Foam::singleComponentSatProp::satL(scalar p) const
{
    notImplemented("bool Foam::singleComponentSatProp::satL(scalar p)");
    return 0.0;
}

// Foam::tmp<Foam::volScalarField>
// Foam::singleComponentSatProp::TSat(const volScalarField& p) const
// {
//     tmp<volScalarField> TSatPtr
//     (
//         new volScalarField
//         (
//             IOobject
//             (
//                 "Tsat",
//                 p.mesh().time().timeName(),
//                 p.mesh(),
//                 IOobject::NO_READ,
//                 IOobject::NO_WRITE
//             ),
//             p.mesh(),
//             dimensionSet(0, 0, 0, 1, 0, 0, 0),
//             calculatedFvPatchScalarField::typeName
//         )
//     );

//     volScalarField& Tsat = TSatPtr.ref();

//     forAll (Tsat.boundaryField(), patchi)
//     {
//         forAll (Tsat.boundaryField()[patchi], i)
//         {
//             Tsat.boundaryFieldRef()[patchi][i] =
//             satT
//             (
//                 p.boundaryField()[patchi][i]
//             );
//         }
//     }

//     forAll (Tsat, celli)
//     {
//         Tsat[celli] = satT(p[celli]);
//     }

//     return TSatPtr;
// }

// Foam::tmp<Foam::volScalarField>
// Foam::singleComponentSatProp::L(const volScalarField& p) const
// {
//     tmp<volScalarField> LPtr
//     (
//         new volScalarField
//         (
//             IOobject
//             (
//                 "L",
//                 p.mesh().time().timeName(),
//                 p.mesh(),
//                 IOobject::NO_READ,
//                 IOobject::NO_WRITE
//             ),
//             p.mesh(),
//             dimEnergy,
//             calculatedFvPatchScalarField::typeName
//         )
//     );

//     volScalarField& LRef = LPtr.ref();

//     forAll (LRef.boundaryField(), patchi)
//     {
//         forAll (LRef.boundaryField()[patchi], i)
//         {
//             LRef.boundaryFieldRef()[patchi][i] =
//             satL
//             (
//                 p.boundaryField()[patchi][i]
//             );
//         }
//     }

//     forAll (LRef, celli)
//     {
//         LRef[celli] = satL(p[celli]);
//     }

//     return LPtr;
// }

// Foam::tmp<Foam::volScalarField>
// Foam::singleComponentSatProp::pSat(const volScalarField& T) const
// {
//     tmp<volScalarField> pSatPtr
//     (
//         new volScalarField
//         (
//             IOobject
//             (
//                 "pSat",
//                 T.mesh().time().timeName(),
//                 T.mesh(),
//                 IOobject::NO_READ,
//                 IOobject::NO_WRITE
//             ),
//             T.mesh(),
//             dimensionSet(0, 2, -2, 0, 0, 0, 0),
//             calculatedFvPatchScalarField::typeName
//         )
//     );

//     volScalarField& pSat = pSatPtr.ref();

//     forAll (pSat.boundaryField(), patchi)
//     {
//         forAll (pSat.boundaryField()[patchi], i)
//         {
//             pSat.boundaryFieldRef()[patchi][i] =
//             satP
//             (
//                 T.boundaryField()[patchi][i]
//             );
//         }
//     }

//     forAll (pSat, celli)
//     {
//         pSat[celli] = satP(T[celli]);
//     }

//     return pSatPtr;
// }

const Foam::dictionary&
Foam::singleComponentSatProp::modelDict() const
{
    return singleComponentSatPropCoeffs_;
}

Foam::dictionary&
Foam::singleComponentSatProp::modelDict()
{
    return singleComponentSatPropCoeffs_;
}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

// bool Foam::singleComponentSatProp::read()
// {
//     return true;
// }

// bool Foam::singleComponentSatProp::correct()
// {

//     volScalarField TLim (max(min(T_, Tmax_), Tmin_));

//     forAll(pSat_, celli)
//     {
//         pSat_[celli] = satP(TLim[celli]);
//     }

//     forAll(pSat_.boundaryField(), patchi)
//     {
//     const fvPatchScalarField& pT = TLim.boundaryField()[patchi];
//     fvPatchScalarField& pPs = pSat_.boundaryFieldRef()[patchi];
//     forAll (pT, i)
//     {
//         pPs[i] = satP(pT[i]);
//     }
//     }

//     return true;
// }


// ************************************************************************* //
