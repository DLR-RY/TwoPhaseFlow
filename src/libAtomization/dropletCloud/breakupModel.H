/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
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
    Foam::breakupModel

Description
    Class for different breakup models for the use with droplet class.

Author
    Dr. Martin Heinrich is copyright owner of the code, based on 
    the breakup models included in the lagrangianSpray library in OpenFOAM.

\*---------------------------------------------------------------------------*/

#ifndef breakupModel_H
#define breakupModel_H

#include "Cloud.H"
#include "dropletCloud.H"
#include "droplet.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class fvMesh;

/*---------------------------------------------------------------------------*\
                           Class dropletCloud Declaration
\*---------------------------------------------------------------------------*/

class breakupModel
{
    // Private data

        const fvMesh& mesh_;

        IOdictionary dict_;

        // Secondary breakup model
        word model_;

        // ETAB model constants
        scalar k1_;
        scalar k2_;
        scalar WeTransition_;
        scalar AWe_;
        scalar TABComega_;
        scalar TABCmu_;
        scalar TABtwoWeCrit_;

        // Reitz Diwakar model constants
        scalar Cbag_;
        scalar Cb_;
        scalar Cstrip_;
        scalar Cs_;

        // Pilch Erman model constants
        scalar B1_;
        scalar B2_;

        // Check for breakup
        void breakup
        (
            const scalar dt,
            const vector g,
            scalar& d,
            scalar& nParticle,
            scalar& y,
            scalar& yDot,
            const scalar rhop,
            const scalar mup,
            const scalar sigma,
            const vector& Up,
            const scalar rhoc,
            const scalar muc,
            const scalar Urmag
        );

public:

    // Constructors

        breakupModel
        (
            const fvMesh&
        );


    // Member Functions

        //- Calulate droplet breakup
        void update
        (
            dropletCloud& cloud,
            droplet::trackingData& td,
            const scalar dt
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
