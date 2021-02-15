/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012-2017 OpenFOAM Foundation
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

#include "solidTabulatedThermo.H"
#include "specie.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class EquationOfState>
inline Foam::solidTabulatedThermo<EquationOfState>::solidTabulatedThermo
(
    const word& name,
    const solidTabulatedThermo& jt
)
:
    EquationOfState(name, jt),
    interpolCpTable_(jt.interpolCpTable_),
    Hf_(jt.Hf_)
{}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class EquationOfState>
inline Foam::solidTabulatedThermo<EquationOfState>::solidTabulatedThermo
(
    const EquationOfState& st,
    const interpolationTable<scalar> cp,
    const scalar hf
)
:
    EquationOfState(st),
    interpolCpTable_(cp),
    Hf_(hf)
{}


template<class EquationOfState>
inline Foam::autoPtr<Foam::solidTabulatedThermo<EquationOfState>>
Foam::solidTabulatedThermo<EquationOfState>::clone() const
{
    return autoPtr<solidTabulatedThermo<EquationOfState>>::New(*this);
}


template<class EquationOfState>
inline Foam::autoPtr<Foam::solidTabulatedThermo<EquationOfState>>
Foam::solidTabulatedThermo<EquationOfState>::New(const dictionary& dict)
{
    return autoPtr<solidTabulatedThermo<EquationOfState>>::New(dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class EquationOfState>
inline Foam::scalar Foam::solidTabulatedThermo<EquationOfState>::limit
(
    const scalar T
) const
{
    return T;
}


template<class EquationOfState>
inline Foam::scalar Foam::solidTabulatedThermo<EquationOfState>::Cp
(
    const scalar p, const scalar T
) const
{
    return interpolCpTable_(T)*this->W();
}


template<class EquationOfState>
inline Foam::scalar Foam::solidTabulatedThermo<EquationOfState>::Ha
(
    const scalar p, const scalar T
) const
{
    return Hs(p, T) + Hc();
}


template<class EquationOfState>
inline Foam::scalar Foam::solidTabulatedThermo<EquationOfState>::Hs
(
    const scalar p, const scalar T
) const
{
    return interpolCpTable_(T)*T*this->W()
            + EquationOfState::H(p, T);
}


template<class EquationOfState>
inline Foam::scalar Foam::solidTabulatedThermo<EquationOfState>::Hc() const
{
    return Hf_;
}


template<class EquationOfState>
inline Foam::scalar Foam::solidTabulatedThermo<EquationOfState>::S
(
    const scalar p, const scalar T
) const
{
    notImplemented
    (
        "scalar solidTabulatedThermo<equationOfState>::"
        "s(const scalar p, const scalar T) const"
    );
    return T;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

// template<class EquationOfState>
// inline void Foam::solidTabulatedThermo<EquationOfState>::operator+=
// (
    // const solidTabulatedThermo<EquationOfState>& ct
// )
// {
    // scalar Y1 = this->Y();

    // EquationOfState::operator+=(ct);

    // if (mag(this->Y()) > SMALL)
    // {
    //     Y1 /= this->Y();
    //     const scalar Y2 = ct.Y()/this->Y();

    //     interpolCpTable_ = Y1*interpolCpTable_ + Y2*ct.interpolCpTable_;
    //     Hf_ = Y1*Hf_ + Y2*ct.Hf_;
    // }
// }


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template<class EquationOfState>
inline Foam::solidTabulatedThermo<EquationOfState> Foam::operator+
(
    const solidTabulatedThermo<EquationOfState>& ct1,
    const solidTabulatedThermo<EquationOfState>& ct2
)
{
    EquationOfState eofs
    (
        static_cast<const EquationOfState&>(ct1)
      + static_cast<const EquationOfState&>(ct2)
    );

    if (mag(eofs.Y()) < SMALL)
    {
        return solidTabulatedThermo<EquationOfState>
        (
            eofs,
            ct1.interpolCpTable_,
            ct1.Hf_
        );
    }
    else
    {
        return solidTabulatedThermo<EquationOfState>
        (
            eofs,
            ct1.Y()/eofs.Y()*ct1.interpolCpTable_
          + ct2.Y()/eofs.Y()*ct2.interpolCpTable_,
            ct1.Y()/eofs.Y()*ct1.Hf_
          + ct2.Y()/eofs.Y()*ct2.Hf_
        );
    }
}


template<class EquationOfState>
inline Foam::solidTabulatedThermo<EquationOfState> Foam::operator*
(
    const scalar s,
    const solidTabulatedThermo<EquationOfState>& ct
)
{
    return solidTabulatedThermo<EquationOfState>
    (
        s*static_cast<const EquationOfState&>(ct),
        ct.interpolCpTable_,
        ct.Hf_
    );
}


template<class EquationOfState>
inline Foam::solidTabulatedThermo<EquationOfState> Foam::operator==
(
    const solidTabulatedThermo<EquationOfState>& ct1,
    const solidTabulatedThermo<EquationOfState>& ct2
)
{
    EquationOfState eofs
    (
        static_cast<const EquationOfState&>(ct1)
     == static_cast<const EquationOfState&>(ct2)
    );

    return solidTabulatedThermo<EquationOfState>
    (
        eofs,
        ct2.Y()/eofs.Y()*ct2.interpolCpTable_
      - ct1.Y()/eofs.Y()*ct1.interpolCpTable_,
        ct2.Y()/eofs.Y()*ct2.Hf_
      - ct1.Y()/eofs.Y()*ct1.Hf_
    );
}


// ************************************************************************* //
