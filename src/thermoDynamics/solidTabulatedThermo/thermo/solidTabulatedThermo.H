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

Class
    Foam::solidTabulatedThermo

Group
    grpSpecieThermo

Description
    tabulated properties thermodynamics package
    templated into the equationOfState.
SourceFiles
    solidTabulatedThermoI.H
    solidTabulatedThermo.C

\*---------------------------------------------------------------------------*/

#ifndef solidTabulatedThermo_H
#define solidTabulatedThermo_H

#include "interpolationTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

template<class EquationOfState> class solidTabulatedThermo;

template<class EquationOfState>
inline solidTabulatedThermo<EquationOfState> operator+
(
    const solidTabulatedThermo<EquationOfState>&,
    const solidTabulatedThermo<EquationOfState>&
);

template<class EquationOfState>
inline solidTabulatedThermo<EquationOfState> operator*
(
    const scalar,
    const solidTabulatedThermo<EquationOfState>&
);


template<class EquationOfState>
inline solidTabulatedThermo<EquationOfState> operator==
(
    const solidTabulatedThermo<EquationOfState>&,
    const solidTabulatedThermo<EquationOfState>&
);


template<class EquationOfState>
Ostream& operator<<
(
    Ostream&,
    const solidTabulatedThermo<EquationOfState>&
);


/*---------------------------------------------------------------------------*\
                         Class solidTabulatedThermo Declaration
\*---------------------------------------------------------------------------*/

template<class EquationOfState>
class solidTabulatedThermo
:
    public EquationOfState
{
    // Private data

        interpolationTable<scalar> interpolCpTable_;
        scalar Hf_;


    // Private Member Functions

        //- Construct from components
        inline solidTabulatedThermo
        (
            const EquationOfState& st,
            const interpolationTable<scalar> interpolCpTable,
            const scalar hf
        );


public:

    // Constructors

        //- Construct from dictionary
        solidTabulatedThermo(const dictionary&);

        //- Construct as a named copy
        inline solidTabulatedThermo
        (
            const word&,
            const solidTabulatedThermo&
        );

         //- Construct and return a clone
        inline autoPtr<solidTabulatedThermo> clone() const;

        //- Selector from dictionary
        inline static autoPtr<solidTabulatedThermo> New(const dictionary& dict);


    // Member Functions

        //- Return the instantiated type name
        static word typeName()
        {
            return "tabulated<" + EquationOfState::typeName() + '>';
        }

        //- Limit the temperature to be in the range Tlow_ to Thigh_
        inline scalar limit(const scalar T) const;


        // Fundamental properties

            //- Heat capacity at constant pressure [J/(kg K)]
            inline scalar Cp(const scalar p, const scalar T) const;

            //- Absolute Enthalpy [J/kg]
            inline scalar Ha(const scalar p, const scalar T) const;

            //- Sensible enthalpy [J/kg]
            inline scalar Hs(const scalar p, const scalar T) const;

            //- Chemical enthalpy [J/kg]
            inline scalar Hc() const;

            //- Entropy [J/(kg K)]
            inline scalar S(const scalar p, const scalar T) const;


    // Member operators

        // inline void operator+=(const solidTabulatedThermo&);


    // Friend operators

        friend solidTabulatedThermo operator+ <EquationOfState>
        (
            const solidTabulatedThermo&,
            const solidTabulatedThermo&
        );

        friend solidTabulatedThermo operator* <EquationOfState>
        (
            const scalar,
            const solidTabulatedThermo&
        );


        friend solidTabulatedThermo operator== <EquationOfState>
        (
            const solidTabulatedThermo&,
            const solidTabulatedThermo&
        );


    // Ostream Operator

        friend Ostream& operator<< <EquationOfState>
        (
            Ostream&,
            const solidTabulatedThermo&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "solidTabulatedThermoI.H"
    #include "solidTabulatedThermo.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
