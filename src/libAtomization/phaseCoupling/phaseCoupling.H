/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Class
    Foam::phaseCoupling

Description
    Class for the coupling between VoF elements and lagrangian droplets
    in the corresponding dropletCloud.

Author
    Dr. Martin Heinrich, all rights reserved

\*---------------------------------------------------------------------------*/

#ifndef phaseCoupling_H
#define phaseCoupling_H

#include "fvCFD.H"
#include "vectorList.H"
#include "dropletCloud.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class phaseCoupling Declaration
\*---------------------------------------------------------------------------*/

class phaseCoupling
{
    // Private data

        //- Reference to mesh and dictionary
        const fvMesh& mesh_;
        const IOdictionary dict_;

        // Reference to velocity field
        const volVectorField& U_;

        // Reference to phase fraction field
        volScalarField& alpha_;

        // Reference to lagrangian droplet cloud from solver
        dropletCloud& cloud_;

	// On/off switch
        bool active_;

        // Limit for alpha which counts towards droplets
        scalar alphaLimit_;

        // Lower and upper limit for diameter to inject droplets
        scalar dMax_;

        // Upper limit for sphericity to inject droplets
        scalar sphMax_;

        // Start time
        scalar startTime_;

        // Interval for phase coupling
        label nInterval_;

        // Marker field to indicate droplet volume ID
        volScalarField vofID_;

        // Field for level set function
        volScalarField psi_;

        // Damping and source functions for droplet removal/creation
        volVectorField source_;
        volScalarField damping_;

        //- Check, if neighbours of cellI belong to connected volI
        void checkConnection
        (
            const label cellI,
            const label iter
        );

        void checkParallelConnection
        (
            const label cellI,
            const label iter
        );
        //- Function to calculate distance to nearest interface
        void levelSetFunction();

public:

    // Constructors

        //- Construct from components
        phaseCoupling
        (
            const volVectorField& U,
            volScalarField& alpha,
            dropletCloud& cloud_
        );


    // Destructor - default

    // Member Functions

        //- Detect droplets and convert to droplets
        void update();

        //- Return source for velocity field
        volVectorField source()
        {
            return source_;
        }

        //- Return damping for velocity field
        volScalarField damping()
        {
            return damping_;
        }


};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
