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
    Foam::dropletCloud

Description
    A cloud of droplets for the VoF-Lagrangian coupling.

Author
    Dr. Martin Heinrich is copyright owner of the code, which is based on 
    the solidParticleCloud class included in the solidParticle library in
    OpenFOAM.

SourceFiles
    dropletCloudI.H
    dropletCloud.C

\*---------------------------------------------------------------------------*/

#ifndef dropletCloud_H
#define dropletCloud_H

#include "Cloud.H"
#include "droplet.H"
#include "collisionModel.H"
#include "breakupModel.H"
#include "IOdictionary.H"
#include "CompactListList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class fvMesh;

/*---------------------------------------------------------------------------*\
                           Class dropletCloud Declaration
\*---------------------------------------------------------------------------*/

class dropletCloud
:
    public Cloud<droplet>
{
    // Private data

        const fvMesh& mesh_;

        const dimensionedVector& g_;
        
        IOdictionary dict_;

        // Droplet properties
        word phaseName_;
        scalar rhop_;
        scalar mup_;
        scalar sigma_;

        //- Face zone IDs for postprocessing
        labelList faceZoneIDs_;

        //- Word list for face zone names
        wordList faceZoneNames_;

        //- Storage for output stream for postprocessing
        List<autoPtr<OFstream> > cloudDataFilesPtr_;
        
        //- Storage for droplet information for postprocessing
        List<List<DynamicList<scalar> > > dataDiameter_;
        List<List<DynamicList<scalar> > > dataNParticle_;
        List<List<DynamicList<vector> > > dataPosition_;
        
        // Auxiliary field for momentum source
        vectorField source_;
        
        // Momentum source for governing equation
        volVectorField momentumSource_;

        //- Class for collision calculation
        collisionModel collision_;

        //- Class for breakup calculation
        breakupModel breakup_;

    // Private Member Functions

        //- Create file for droplet postprocessing
        void createDropletDataFiles();

        //- Give a short summary of droplet cloud
        void info();

        //- Calculate mean diameter
        Foam::scalar Dij(const label i, const label j);

        //- Calculate maximum diameter
        Foam::scalar Dmax();


public:

    // Constructors

        //- Construct given mesh
        dropletCloud
        (
            const fvMesh&,
            const dimensionedVector&, 
            const word& cloudName = "dropletCloud",
            bool readFields = true
        );


    // Member Functions

        // Access

            inline const fvMesh& mesh() const;

            inline scalar rhop() const;
            inline scalar e() const;
            inline scalar mup() const;
            inline scalar sigma() const;
            
            // Field to store momentum source
            inline vectorField& source();
            
            // Return momentum source to carrier phase
            inline const volVectorField& momentumSource();
            
            // Return faceZoneIDs for hitFace in droplet
            inline const labelList& faceZoneIDs();
            
            // Return processor data for droplet postprocessing
            inline List<List<DynamicList<scalar> > >& dataDiameter();
            inline List<List<DynamicList<vector> > >& dataPosition();
            inline List<List<DynamicList<scalar> > >& dataNParticle();


        // Edit

            //- Move the droplets under the influence of the given
            //  gravitational acceleration
            void move();
  
            //- Inject droplets at a given position with a given
            //  velocity and diameter
            void inject
            (
                vector position,
                const scalar diameter,
                const vector velocity
            );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "dropletCloudI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
