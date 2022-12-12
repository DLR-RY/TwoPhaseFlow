/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
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

#include "makeThermo.H"
#include "makeReactionThermo.H"

#include "rhoThermo.H"
#include "rhoReactionThermo.H"


#include "specie.H"
#include "perfectGas.H"
#include "PengRobinsonGas.H"

#include "hConstThermo.H"
#include "sutherlandTransport.H"

#include "sensibleInternalEnergy.H"
#include "thermo.H"

#include "constTransport.H"
#include "hPolynomialThermo.H"
#include "polynomialTransport.H"

// #include "hTabulatedThermo.H"
// #include "tabulatedTransport.H"

#include "heRhoThermo.H"
#include "pureMixture.H"
#include "homogeneousMixture.H"
#include "multiComponentMixture.H"

// #include "thermoPhysicsTypes.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// added polyThermo
makeThermo
(
    rhoThermo,
    heRhoThermo,
    pureMixture,
    polynomialTransport,
    sensibleInternalEnergy,
    hPolynomialThermo,
    perfectGas,
    specie
);

makeThermos
(
    rhoThermo,
    heRhoThermo,
    pureMixture,
    polynomialTransport,
    sensibleInternalEnergy,
    hPolynomialThermo,
    PengRobinsonGas,
    specie
);

// added tabThermo
// makeThermo
// (
//     rhoThermo,
//     heRhoThermo,
//     pureMixture,
//     tabulatedTransport,
//     sensibleInternalEnergy,
//     hTabulatedThermo,
//     perfectGas,
//     specie
// );

// makeThermo
// (
//     rhoThermo,
//     heRhoThermo,
//     pureMixture,
//     tabulatedTransport,
//     sensibleInternalEnergy,
//     hTabulatedThermo,
//     PengRobinsonGas,
//     specie
// );


// thermo physics types based on sensibleInternalEnergy
// typedef
// tabulatedTransport
// <
//     species::thermo
//     <
//         hTabulatedThermo
//         <
//             perfectGas<specie>
//         >,
//         sensibleInternalEnergy
//     >
// > tabTabEThermoPhysics;

// typedef
// tabulatedTransport
// <
//     species::thermo
//     <
//         hTabulatedThermo
//         <
//             PengRobinsonGas<specie>
//         >,
//         sensibleInternalEnergy
//     >
// > tabTabPengRobPhysics;

typedef
polynomialTransport
<
    species::thermo
    <
        hPolynomialThermo
        <
            perfectGas<specie>,
            8
        >,
        sensibleInternalEnergy
    >,
    8
> polyPolyEThermoPhysics;

makeThermoPhysicsReactionThermos
(
    rhoThermo,
    rhoReactionThermo,
    heRhoThermo,
    multiComponentMixture,
    polyPolyEThermoPhysics
);

// makeThermoPhysicsReactionThermos
// (
//     rhoThermo,
//     rhoReactionThermo,
//     heRhoThermo,
//     multiComponentMixture,
//     tabTabEThermoPhysics
// );


// makeThermoPhysicsReactionThermos
// (
//     rhoThermo,
//     rhoReactionThermo,
//     heRhoThermo,
//     multiComponentMixture,
//     tabTabPengRobPhysics
// );

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
