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

#include "hTabulatedThermo.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class EquationOfState>
inline Foam::hTabulatedThermo<EquationOfState>::hTabulatedThermo
(
    const EquationOfState& pt,
    const scalar Hf,
    const scalar Sf,
    const Pair<scalarField>& CpTable_
)
:
    EquationOfState(pt),
    Hf_(Hf),
    Sf_(Sf),
    CpTable_(CpTable_)
{}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class EquationOfState>
inline Foam::hTabulatedThermo<EquationOfState>::hTabulatedThermo
(
    const word& name,
    const hTabulatedThermo& pt
)
:
    EquationOfState(name, pt),
    Hf_(pt.Hf_),
    Sf_(pt.Sf_),
    CpTable_(pt.CpTable_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class EquationOfState>
inline Foam::scalar Foam::hTabulatedThermo<EquationOfState>::limit
(
    const scalar T
) const
{
    return T;
}


template<class EquationOfState>
inline Foam::scalar Foam::hTabulatedThermo<EquationOfState>::Cp
(
    const scalar p, const scalar T
) const
{

    return interpolateXY(T,CpTable_.first(),CpTable_.second())
            + EquationOfState::Cp(p, T);
}


template<class EquationOfState>
inline Foam::scalar Foam::hTabulatedThermo<EquationOfState>::Ha
(
    const scalar p, const scalar T
) const
{
    return interpolateXY(T,CpTable_.first(),CpTable_.second())*T
             + EquationOfState::H(p, T);
}


template<class EquationOfState>
inline Foam::scalar Foam::hTabulatedThermo<EquationOfState>::Hs
(
    const scalar p, const scalar T
) const
{
    return Ha(p, T) - Hc();
}


template<class EquationOfState>
inline Foam::scalar Foam::hTabulatedThermo<EquationOfState>::Hc()
const
{
    return Hf_;
}


template<class EquationOfState>
inline Foam::scalar Foam::hTabulatedThermo<EquationOfState>::S
(
    const scalar p,
    const scalar T
) const
{
    notImplemented("S")
    return 0; //sCoeffs_.value(T) + EquationOfState::S(p, T);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class EquationOfState>
inline void Foam::hTabulatedThermo<EquationOfState>::operator=
(
    const hTabulatedThermo<EquationOfState>& pt
)
{
    EquationOfState::operator=(pt);

    Hf_ = pt.Hf_;
    Sf_ = pt.Sf_;
    CpTable_ = pt.CpTable_;
}


template<class EquationOfState>
inline void Foam::hTabulatedThermo<EquationOfState>::operator+=
(
    const hTabulatedThermo<EquationOfState>& pt
)
{
    scalar Y1 = this->Y();

    EquationOfState::operator+=(pt);

    if (mag(this->Y()) > SMALL)
    {
        Y1 /= this->Y();
        const scalar Y2 = pt.Y()/this->Y();

        Hf_ = Y1*Hf_ + Y2*pt.Hf_;
        Sf_ = Y1*Sf_ + Y2*pt.Sf_;
        if
        (
            CpTable_.first().size() == pt.CpTable_.first().size()
        &&  CpTable_.first().first() == pt.CpTable_.first().first()
        &&  CpTable_.first().last() == pt.CpTable_.first().last()
        )
        {
            CpTable_.second() = Y1*CpTable_.second() + Y2*pt.CpTable_.second();
        }
        else
        {
            FatalErrorInFunction
                << "sample rate,Tmin and Tmax have to be identical"
                << exit(FatalError);
        }

        // CpTable_ = pt.CpTable_;
        // notImplemented("cannot weight tables linearly")

    }
}


template<class EquationOfState>
inline void Foam::hTabulatedThermo<EquationOfState>::operator*=
(
    const scalar s
)
{
    EquationOfState::operator*=(s);
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template<class EquationOfState>
inline Foam::hTabulatedThermo<EquationOfState> Foam::operator+
(
    const hTabulatedThermo<EquationOfState>& pt1,
    const hTabulatedThermo<EquationOfState>& pt2
)
{
    EquationOfState eofs = pt1;
    eofs += pt2;

    if (mag(eofs.Y()) < SMALL)
    {
        return hTabulatedThermo<EquationOfState>
        (
            eofs,
            pt1.Hf_,
            pt1.Sf_,
            pt1.CpTable_
        );
    }
    {
        const scalar Y1 = pt1.Y()/eofs.Y();
        const scalar Y2 = pt2.Y()/eofs.Y();

        return hTabulatedThermo<EquationOfState>
        (
            eofs,
            Y1*pt1.Hf_ + Y2*pt2.Hf_,
            Y1*pt1.Sf_ + Y2*pt2.Sf_,
            Y1*pt1.CpTable_ + Y2*pt2.CpTable_ // makes no sense
        );
    }
}


template<class EquationOfState>
inline Foam::hTabulatedThermo<EquationOfState> Foam::operator*
(
    const scalar s,
    const hTabulatedThermo<EquationOfState>& pt
)
{
    return hTabulatedThermo<EquationOfState>
    (
        s*static_cast<const EquationOfState&>(pt),
        pt.Hf_,
        pt.Sf_,
        pt.CpTable_
    );
}


template<class EquationOfState>
inline Foam::hTabulatedThermo<EquationOfState> Foam::operator==
(
    const hTabulatedThermo<EquationOfState>& pt1,
    const hTabulatedThermo<EquationOfState>& pt2
)
{
    EquationOfState eofs
    (
        static_cast<const EquationOfState&>(pt1)
     == static_cast<const EquationOfState&>(pt2)
    );

    const scalar Y1 = pt1.Y()/eofs.Y();
    const scalar Y2 = pt2.Y()/eofs.Y();

    return hTabulatedThermo<EquationOfState>
    (
        eofs,
        Y2*pt2.Hf_       - Y1*pt1.Hf_,
        Y2*pt2.Sf_       - Y1*pt1.Sf_,
        Y1*pt1.CpTable_ - Y2*pt2.CpTable_ // makes no sense
    );
}


// ************************************************************************* //
