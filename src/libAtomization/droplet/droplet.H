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
    Foam::droplet

Description
    Spherical droplet class with optional two-way coupling with the
    continuous phase, droplet-droplet collision, and breakup..

Author
    Dr. Martin Heinrich is copyright owner of the additional code, modified
    on the basis of the solidParticle class which is part of OpenFOAM.

SourceFiles
    dropletI.H
    droplet.C
    dropletIO.C

\*---------------------------------------------------------------------------*/

#ifndef droplet_H
#define droplet_H

#include "particle.H"
#include "IOstream.H"
#include "autoPtr.H"
#include "interpolationCellPoint.H"
#include "contiguous.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class dropletCloud;


// Forward declaration of friend functions and operators

class droplet;

Ostream& operator<<(Ostream&, const droplet&);


/*---------------------------------------------------------------------------*\
                           Class droplet Declaration
\*---------------------------------------------------------------------------*/

class droplet
:
    public particle
{
    // Private data

        //- Diameter
        scalar d_;

        //- Velocity of particle
        vector U_;

        //- Number of particles in parcel
        scalar nParticle_;

        //- Spherical deviation
	scalar y_;

	//- Rate of change of spherical deviation
	scalar yDot_;

public:

    friend class Cloud<droplet>;

    //- Class used to pass tracking data to the trackToFace function
    class trackingData
    :
        public particle::trackingData
    {
        // Interpolators for continuous phase fields

            const interpolationCellPoint<scalar>& rhoInterp_;
            const interpolationCellPoint<vector>& UInterp_;
            const interpolationCellPoint<scalar>& nuInterp_;

        //- Local gravitational or other body-force acceleration
        const vector& g_;


    public:

        // Constructors

            inline trackingData
            (
                const dropletCloud& spc,
                const interpolationCellPoint<scalar>& rhoInterp,
                const interpolationCellPoint<vector>& UInterp,
                const interpolationCellPoint<scalar>& nuInterp,
                const vector& g
            );


        // Member functions

            inline const interpolationCellPoint<scalar>& rhoInterp() const;

            inline const interpolationCellPoint<vector>& UInterp() const;

            inline const interpolationCellPoint<scalar>& nuInterp() const;

            inline const vector& g() const;
    };


    // Static data members

        //- Size in bytes of the fields
        static const std::size_t sizeofFields;


    // Constructors

        //- Construct from components
        inline droplet
        (
            const polyMesh& mesh,
            const vector& position,
            const label celli,
            const scalar d,
            const vector& U
        );

        //- Construct from Istream
        droplet
        (
            const polyMesh& mesh,
            Istream& is,
            bool readFields = true,
            bool newFormat = true
        );

        //- Construct and return a clone
        virtual autoPtr<particle> clone() const
        {
            return autoPtr<particle>(new droplet(*this));
        }

        //- Factory class to read-construct particles used for
        //  parallel transfer
        class iNew
        {
            const polyMesh& mesh_;

        public:

            iNew(const polyMesh& mesh)
            :
                mesh_(mesh)
            {}

            autoPtr<droplet> operator()(Istream& is) const
            {
                return autoPtr<droplet>
                (
                    new droplet(mesh_, is, true)
                );
            }
        };


    // Member Functions

        // Access

            //- Return diameter
            inline scalar& d();

            //- Return velocity
            inline vector& U();

            //- Return number of particles in parcel
            inline scalar& nParticle();

            //- Return spherical deviation
            inline scalar& y();

            //- Return change of rate of spherical deviation
            inline scalar& yDot();


        // Tracking

            //- Move
            bool move(dropletCloud&, trackingData&, const scalar);


        // Patch interactions

            //- Overridable function to handle the particle hitting a patch
            //  Executed before other patch-hitting functions
            bool hitPatch(dropletCloud& cloud, trackingData& td);

            //- Overridable function to handle the particle hitting a
            //  processorPatch
            void hitProcessorPatch(dropletCloud& cloud, trackingData& td);

            //- Overridable function to handle the particle hitting a wallPatch
            void hitWallPatch(dropletCloud& cloud, trackingData& td);

            //- Transform the physical properties of the particle
            //  according to the given transformation tensor
            virtual void transformProperties(const tensor& T);

            //- Transform the physical properties of the particle
            //  according to the given separation vector
            virtual void transformProperties(const vector& separation);


    // I-O

        static void readFields(Cloud<droplet>& c);

        static void writeFields(const Cloud<droplet>& c);


    // Ostream Operator

        friend Ostream& operator<<(Ostream&, const droplet&);
};

//- Contiguous data for solidParticle
template<> struct is_contiguous<droplet> : std::true_type {};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "dropletI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
