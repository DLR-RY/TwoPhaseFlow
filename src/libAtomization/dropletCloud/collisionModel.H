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
    Foam::collisionModel

Description
    Class for the collision model for the use with droplet class.

Author
    Dr. Martin Heinrich is copyright owner of the code, based on 
    the collision models included in the lagrangianSpray library in OpenFOAM.


\*---------------------------------------------------------------------------*/

#ifndef collisionModel_H
#define collisionModel_H

#include "Cloud.H"
#include "dropletCloud.H"
#include "droplet.H"
#include "Random.H"
#include "clock.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class fvMesh;

/*---------------------------------------------------------------------------*\
                           Class dropletCloud Declaration
\*---------------------------------------------------------------------------*/

class collisionModel
{
    // Private data

        const fvMesh& mesh_;

	IOdictionary dict_;

	// Secondary breakup active
	bool active_;

	//- Random number generator
        Random ranGen_;

        //- Model constants
        scalar cTime_;
        scalar cSpace_;

        //- Droplet density
        scalar rhop_;

        //- Droplet surface tension
        scalar sigma_;

    // Private Member Functions

        //- Collide parcels and return true if mass has changed
        bool collideDroplets
        (
            const scalar dt,
            droplet& p1,
            droplet& p2,
            scalar& m1,
            scalar& m2
        );

        // 1 is the larger drop and 2 is the smaller
        bool collideSorted
        (
            const scalar dt,
            droplet& p1,
            droplet& p2,
            scalar& m1,
            scalar& m2
        );
           

public:

    // Constructors

        collisionModel
        (
            const fvMesh&
        );


    // Member Functions

        //- Calulate droplet-droplet collision and coalescence
        void update
        (
            dropletCloud& cloud,
            const scalar dt
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
