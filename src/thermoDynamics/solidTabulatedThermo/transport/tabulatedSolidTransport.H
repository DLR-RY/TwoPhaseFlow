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
    Foam::tabulatedSolidTransport

Description
    Constant properties Transport package.
    Templated into a given thermodynamics package (needed for thermal
    conductivity).

SourceFiles
    tabulatedSolidTransportI.H
    tabulatedSolidTransport.C

\*---------------------------------------------------------------------------*/

#ifndef tabulatedSolidTransport_H
#define tabulatedSolidTransport_H

#include "vector.H"
#include "interpolationTable.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<class Thermo> class tabulatedSolidTransport;

template<class Thermo>
inline tabulatedSolidTransport<Thermo> operator*
(
    const scalar,
    const tabulatedSolidTransport<Thermo>&
);

template<class Thermo>
Ostream& operator<<
(
    Ostream&,
    const tabulatedSolidTransport<Thermo>&
);


/*---------------------------------------------------------------------------*\
                   Class tabulatedSolidTransport Declaration
\*---------------------------------------------------------------------------*/

template<class Thermo>
class tabulatedSolidTransport
:
    public Thermo
{
    // Private data

        //- Constant isotropic thermal conductivity
//        scalar kappa_;
        interpolationTable<scalar> interpolCpTable_;


    // Private Member Functions

        //- Construct from components
        inline tabulatedSolidTransport(const Thermo& t, const interpolationTable<scalar> interpolCpTable);


public:

    // Constructors

        //- Construct as named copy
        inline tabulatedSolidTransport
        (
            const word&,
            const tabulatedSolidTransport&
        );

        //- Construct from dictionary
        tabulatedSolidTransport(const dictionary& dict);

        // Selector from dictionary
        inline static autoPtr<tabulatedSolidTransport> New
        (
            const dictionary& dict
        );


    // Member functions

        //- Return the instantiated type name
        static word typeName()
        {
            return "isoTabulated<" + Thermo::typeName() + '>';
        }

        //- Is the thermal conductivity isotropic
        static const bool isotropic = true;

        //- Isotropic thermal conductivity [W/mK]
        inline scalar kappa(const scalar p, const scalar T) const;

        //- Un-isotropic thermal conductivity [W/mK]
        inline vector Kappa(const scalar p, const scalar T) const;

         //- Dynamic viscosity [kg/ms]
        inline scalar mu(const scalar p, const scalar T) const;

        //- Thermal diffusivity of enthalpy [kg/ms]
        inline scalar alphah(const scalar p, const scalar T) const;


        //- Write to Ostream
        void write(Ostream& os) const;


    // Member operators

        inline void operator=(const tabulatedSolidTransport&);
        // inline void operator+=(const tabulatedSolidTransport&);


    // Friend operators

        friend tabulatedSolidTransport operator* <Thermo>
        (
            const scalar,
            const tabulatedSolidTransport&
        );


    // Ostream Operator

        friend Ostream& operator<< <Thermo>
        (
            Ostream&,
            const tabulatedSolidTransport&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "tabulatedSolidTransportI.H"

#ifdef NoRepository
    #include "tabulatedSolidTransport.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
