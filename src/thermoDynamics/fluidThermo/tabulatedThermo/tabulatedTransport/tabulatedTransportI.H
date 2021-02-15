/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2017 OpenFOAM Foundation
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

#include "specie.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Thermo>
inline Foam::tabulatedTransport<Thermo>::tabulatedTransport
(
    const Thermo& t,
    const Pair<scalarField>& MuTable,
    const Pair<scalarField>& KappaTable
)
:
    Thermo(t),
    MuTable_(MuTable),
    KappaTable_(KappaTable)
{}


template<class Thermo>
inline Foam::tabulatedTransport<Thermo>::tabulatedTransport
(
    const word& name,
    const tabulatedTransport& pt
)
:
    Thermo(name, pt),
    MuTable_(pt.MuTable_),
    KappaTable_(pt.KappaTable_)
{}


template<class Thermo>
inline Foam::autoPtr<Foam::tabulatedTransport<Thermo>>
Foam::tabulatedTransport<Thermo>::clone() const
{
    return autoPtr<tabulatedTransport<Thermo>>::New(*this);
}


template<class Thermo>
inline Foam::autoPtr<Foam::tabulatedTransport<Thermo>>
Foam::tabulatedTransport<Thermo>::New(const dictionary& dict)
{
    return autoPtr<tabulatedTransport<Thermo>>::New(dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Thermo>
inline Foam::scalar Foam::tabulatedTransport<Thermo>::mu
(
    const scalar p,
    const scalar T
) const
{
    return interpolateXY(T,MuTable_.first(),MuTable_.second());
}

template<class Thermo>
inline Foam::scalar Foam::tabulatedTransport<Thermo>::kappa
(
    const scalar p,
    const scalar T
) const
{
    return interpolateXY(T,KappaTable_.first(),KappaTable_.second());
}


template<class Thermo>
inline Foam::scalar Foam::tabulatedTransport<Thermo>::alphah
(
    const scalar p, const scalar T
) const
{
    return kappa(p, T)/this->Cp(p, T);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Thermo>
inline void Foam::tabulatedTransport<Thermo>::operator=
(
    const tabulatedTransport<Thermo>& pt
)
{
    Thermo::operator=(pt);

    MuTable_ = pt.MuTable_;
    KappaTable_ = pt.KappaTable_;
}


template<class Thermo>
inline void Foam::tabulatedTransport<Thermo>::operator+=
(
    const tabulatedTransport<Thermo>& pt
)
{
    scalar Y1 = this->Y();

    Thermo::operator+=(pt);

    if (mag(this->Y()) > SMALL)
    {
        Y1 /= this->Y();
        scalar Y2 = pt.Y()/this->Y();
        // MuTable and KappaTable are sampled identical
        // so we donot have to check both
        if
        (
            MuTable_.first().size() == pt.MuTable_.first().size()
        &&  MuTable_.first().first() == pt.MuTable_.first().first()
        &&  MuTable_.first().last() == pt.MuTable_.first().last()
        )
        {
            MuTable_.second() = Y1*MuTable_.second() + Y2*pt.MuTable_.second();
            KappaTable_.second() =
                Y1*KappaTable_.second() + Y2*pt.KappaTable_.second();
        }
        else
        {
            FatalErrorInFunction
                << "nPoints,Tmin and Tmax have to be identical"
                << exit(FatalError);
        }

        //MuTable_ = Y1*MuTable_ + Y2*pt.MuTable_;
        //KappaTable_ = Y1*KappaTable_ + Y2*pt.KappaTable_;
    }
}


template<class Thermo>
inline void Foam::tabulatedTransport<Thermo>::operator*=
(
    const scalar s
)
{
    Thermo::operator*=(s);
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template<class Thermo>
inline Foam::tabulatedTransport<Thermo> Foam::operator+
(
    const tabulatedTransport<Thermo>& pt1,
    const tabulatedTransport<Thermo>& pt2
)
{
    Thermo t
    (
        static_cast<const Thermo&>(pt1) + static_cast<const Thermo&>(pt2)
    );

    if (mag(t.Y()) < SMALL)
    {
        return tabulatedTransport<Thermo>
        (
            t,
            0,
            pt1.MuTable_,
            pt1.KappaTable_
        );
    }
    else
    {
        scalar Y1 = pt1.Y()/t.Y();
        scalar Y2 = pt2.Y()/t.Y();

        return tabulatedTransport<Thermo>
        (
            t,
            Y1*pt1.MuTable_ + Y2*pt2.MuTable_,
            Y1*pt1.KappaTable_ + Y2*pt2.KappaTable_
        );
    }
}


template<class Thermo>
inline Foam::tabulatedTransport<Thermo> Foam::operator*
(
    const scalar s,
    const tabulatedTransport<Thermo>& pt
)
{
    return tabulatedTransport<Thermo>
    (
        s*static_cast<const Thermo&>(pt),
        pt.MuTable_,
        pt.KappaTable_
    );
}


// ************************************************************************* //
