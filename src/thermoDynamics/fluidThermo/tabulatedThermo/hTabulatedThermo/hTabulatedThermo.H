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

Class
    Foam::hTabulatedThermo

Group
    grpSpecieThermo

Description
    Thermodynamics package templated on the equation of state, using polynomial
    functions for \c cp, \c h and \c s.

    Polynomials for \c h and \c s derived from \c cp.

Usage

    \table
        Property     | Description
        Hf           | Heat of formation
        Sf           | Standard entropy
        CpFile       | Specific heat at constant pressure table
    \endtable

    Example of the specification of the thermodynamic properties:
    \verbatim
    thermodynamics
    {
        Hf              0;
        Sf              0;
        CpFile          "cp.dat"
    }
    \endverbatim

    cat cp.dat:
    2
    (
        (1 2300)
        (100 2300)
    )

Note
    - Heat of formation is inputted in [J/kg], but internally uses [J/kmol]
    - Standard entropy is inputted in [J/kg/K], but internally uses [J/kmol/K]
    - Specific heat at constant pressure polynomial coefficients evaluate to an
      expression in [J/(kg.K)].

SourceFiles
    hTabulatedThermoI.H
    hTabulatedThermo.C

See also
    Foam::Polynomial

\*---------------------------------------------------------------------------*/

#ifndef hTabulatedThermo_H
#define hTabulatedThermo_H

#include "scalar.H"
#include "interpolationTable.H"
#include "interpolateXY.H"
//interpolation2DTable<scalar> interpolCpTable_;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

template<class EquationOfState>
class hTabulatedThermo;

template<class EquationOfState>
inline hTabulatedThermo<EquationOfState> operator+
(
    const hTabulatedThermo<EquationOfState>&,
    const hTabulatedThermo<EquationOfState>&
);

template<class EquationOfState>
inline hTabulatedThermo<EquationOfState> operator*
(
    const scalar,
    const hTabulatedThermo<EquationOfState>&
);

template<class EquationOfState>
inline hTabulatedThermo<EquationOfState> operator==
(
    const hTabulatedThermo<EquationOfState>&,
    const hTabulatedThermo<EquationOfState>&
);

template<class EquationOfState>
Ostream& operator<<
(
    Ostream&,
    const hTabulatedThermo<EquationOfState>&
);


/*---------------------------------------------------------------------------*\
                      Class hTabulatedThermo Declaration
\*---------------------------------------------------------------------------*/

template<class EquationOfState>
class hTabulatedThermo
:
    public EquationOfState
{
    // Private data

        //- Heat of formation
        scalar Hf_;

        //- Standard entropy
        scalar Sf_;

        //- data table
        Pair<scalarField> CpTable_;



    // Private Member Functions

        void readTable
        (
            tableReader<scalar>& reader,
            List<Tuple2<scalar, scalar>>& data,
            fileName file
        );

        void check
        (
            List<Tuple2<scalar, scalar>>& data
        ) const;

        //- Coeffs name. Eg, "CpCoeffs<10>"
        // inline static word coeffsName(const char* name)
        // {
        //     return word(name + ("Coeffs<" + std::to_string(PolySize) + '>'));
        // }

        //- Construct from components
        inline hTabulatedThermo
        (
            const EquationOfState& pt,
            const scalar Hf,
            const scalar Sf,
            const Pair<scalarField>& CpTable_
        );


public:

    // Constructors

        //- Construct from dictionary
        hTabulatedThermo(const dictionary& dict);

        //- Construct as a named copy
        inline hTabulatedThermo(const word&, const hTabulatedThermo&);


    // Member Functions

        //- Return the instantiated type name
        static word typeName()
        {
            return "hTabulated<" + EquationOfState::typeName() + '>';
        }

        //- Limit the temperature to be in the range Tlow_ to Thigh_
        inline scalar limit(const scalar) const;


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


        // I-O

            //- Write to Ostream
            void write(Ostream& os) const;


    // Member operators

        inline void operator=(const hTabulatedThermo&);
        inline void operator+=(const hTabulatedThermo&);
        inline void operator*=(const scalar);


    // Friend operators

        friend hTabulatedThermo operator+ <EquationOfState>
        (
            const hTabulatedThermo&,
            const hTabulatedThermo&
        );

        friend hTabulatedThermo operator* <EquationOfState>
        (
            const scalar,
            const hTabulatedThermo&
        );

        friend hTabulatedThermo operator== <EquationOfState>
        (
            const hTabulatedThermo&,
            const hTabulatedThermo&
        );


    // Ostream Operator

        friend Ostream& operator<< <EquationOfState>
        (
            Ostream&,
            const hTabulatedThermo&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "hTabulatedThermoI.H"

#ifdef NoRepository
    #include "hTabulatedThermo.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
