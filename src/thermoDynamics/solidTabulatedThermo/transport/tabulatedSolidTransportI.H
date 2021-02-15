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

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class thermo>
inline Foam::tabulatedSolidTransport<thermo>::tabulatedSolidTransport
(
    const thermo& t,
    const interpolationTable<scalar> cp
)
:
    thermo(t),
    interpolCpTable_(cp)
{}


template<class thermo>
inline Foam::tabulatedSolidTransport<thermo>::tabulatedSolidTransport
(
    const word& name,
    const tabulatedSolidTransport& ct
)
:
    thermo(name, ct),
    interpolCpTable_(ct.interpolCpTable_)
{}


template<class Thermo>
inline Foam::autoPtr<Foam::tabulatedSolidTransport<Thermo>>
Foam::tabulatedSolidTransport<Thermo>::New
(
    const dictionary& dict
)
{
    return autoPtr<tabulatedSolidTransport<Thermo>>::New(dict);
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class thermo>
inline Foam::scalar Foam::tabulatedSolidTransport<thermo>::
kappa(const scalar p, const scalar T) const
{
    return interpolCpTable_(T);
}

template<class thermo>
inline Foam::vector Foam::tabulatedSolidTransport<thermo>::
Kappa(const scalar p, const scalar T) const
{
    return vector(interpolCpTable_(T), interpolCpTable_(T), interpolCpTable_(T));
}


template<class thermo>
inline Foam::scalar Foam::tabulatedSolidTransport<thermo>::
mu(const scalar p, const scalar T) const
{
    NotImplemented;
    return scalar(0);
}


template<class thermo>
inline Foam::scalar Foam::tabulatedSolidTransport<thermo>::
alphah(const scalar p, const scalar T) const
{
    return interpolCpTable_(T)/this->Cp(p, T);
}

// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class thermo>
inline void Foam::tabulatedSolidTransport<thermo>::operator=
(
    const tabulatedSolidTransport<thermo>& ct
)
{
    thermo::operator=(ct);
    interpolCpTable_ = ct.interpolCpTable_;
}


// template<class thermo>
// inline void Foam::tabulatedSolidTransport<thermo>::operator+=
// (
//     const tabulatedSolidTransport<thermo>& ct
// )
// {
//     scalar Y1 = this->Y();
//     thermo::operator+=(ct);

//     Y1 /= this->Y();
//     scalar Y2 = ct.Y()/this->Y();

//     interpolCpTable_ = Y1*interpolCpTable_ + Y2*ct.interpolCpTable_;
// }


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //


template<class thermo>
inline Foam::tabulatedSolidTransport<thermo> Foam::operator*
(
    const scalar s,
    const tabulatedSolidTransport<thermo>& ct
)
{
    return tabulatedSolidTransport<thermo>
    (
        s*static_cast<const thermo&>(ct),
        ct.interpolCpTable_
    );
}


// ************************************************************************* //
