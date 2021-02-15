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
    Foam::tabulatedTransport

Group
    grpSpecieTransport

Description
    Transport package using tables for \c mu and \c kappa.

Usage

    \table
        Property        | Description
        MuFile          | Dynamic viscosity table
        kappaFile       | Thermal conductivity table
    \endtable

    Example of the specification of the transport properties:
    \verbatim
    transport
    {
        MuFile     "mu.dat"
        KappaFile  "kappa.dat"
    }
    \endverbatim

    cat mu.dat:
    2
    (
        (1 2300)
        (100 2300)
    )

Note
    - Dynamic viscosity polynomial coefficients evaluate to an expression in
      [Pa.s], but internally uses [Pa.s/kmol].
    - Thermal conductivity polynomial coefficients evaluate to an expression in
      [W/m/K], but internally uses [W/m/K/kmol].

SourceFiles
    tabulatedTransportI.H
    tabulatedTransport.C

See also
    Foam::Polynomial

\*---------------------------------------------------------------------------*/

#ifndef tabulatedTransport_H
#define tabulatedTransport_H

#include "interpolateXY.H"
#include "interpolationTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

template<class Thermo> class tabulatedTransport;

template<class Thermo>
inline tabulatedTransport<Thermo> operator+
(
    const tabulatedTransport<Thermo>&,
    const tabulatedTransport<Thermo>&
);

template<class Thermo>
inline tabulatedTransport<Thermo> operator*
(
    const scalar,
    const tabulatedTransport<Thermo>&
);

template<class Thermo>
Ostream& operator<<
(
    Ostream&,
    const tabulatedTransport<Thermo>&
);


/*---------------------------------------------------------------------------*\
                     Class tabulatedTransport Declaration
\*---------------------------------------------------------------------------*/

template<class Thermo>
class tabulatedTransport
:
    public Thermo
{
    // Private data

        //- Dynamic viscosity polynomial coefficients
        // interpolationTable<scalar> interpolMuTable_;
        Pair<scalarField> MuTable_;

        //- Thermal conductivity polynomial coefficients
        // interpolationTable<scalar> interpolKappaTable_;
        Pair<scalarField> KappaTable_;


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

        //- Coeffs name. Eg, "muCoeffs<10>"
        // inline static word coeffsName(const char* name)
        // {
        //     return word(name + ("Coeffs<" + std::to_string(PolySize) + '>'));
        // }

        //- Construct from components
        inline tabulatedTransport
        (
            const Thermo& t,
            const Pair<scalarField>& MuTable,
            const Pair<scalarField>& KappaTable
        );


public:

    // Constructors

        //- Construct as named copy
        inline tabulatedTransport(const word&, const tabulatedTransport&);

        //- Construct from dictionary
        tabulatedTransport(const dictionary& dict);

        //- Construct and return a clone
        inline autoPtr<tabulatedTransport> clone() const;

        // Selector from dictionary
        inline static autoPtr<tabulatedTransport> New(const dictionary& dict);


    // Member functions

        //- Return the instantiated type name
        static word typeName()
        {
            return "tabulated<" + Thermo::typeName() + '>';
        }

        //- Dynamic viscosity [kg/ms]
        inline scalar mu(const scalar p, const scalar T) const;

        //- Thermal conductivity [W/mK]
        inline scalar kappa(const scalar p, const scalar T) const;

        //- Thermal diffusivity of enthalpy [kg/ms]
        inline scalar alphah(const scalar p, const scalar T) const;

        // Species diffusivity
        //inline scalar D(const scalar p, const scalar T) const;

        //- Write to Ostream
        void write(Ostream& os) const;


    // Member operators

        inline void operator=(const tabulatedTransport&);

        inline void operator+=(const tabulatedTransport&);

        inline void operator*=(const scalar);


    // Friend operators

        friend tabulatedTransport operator+ <Thermo>
        (
            const tabulatedTransport&,
            const tabulatedTransport&
        );

        friend tabulatedTransport operator* <Thermo>
        (
            const scalar,
            const tabulatedTransport&
        );


    // Ostream Operator

        friend Ostream& operator<< <Thermo>
        (
            Ostream&,
            const tabulatedTransport&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "tabulatedTransportI.H"

#ifdef NoRepository
    #include "tabulatedTransport.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
