/*---------------------------------------------------------------------------*\
            Copyright (c) 2017-2019, German Aerospace Center (DLR)
-------------------------------------------------------------------------------
License
    This file is part of the VoFLibrary source code library, which is an
    unofficial extension to OpenFOAM.

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
    Foam::sampledInterface

Description

Author
    Henning Scheufler, DLR, all rights reserved.

SourceFiles
    sampledInterface.C

\*---------------------------------------------------------------------------*/

#ifndef sampledInterface_H
#define sampledInterface_H

#include "interface.H"
#include "sampledSurface.H"
#include "ZoneIDs.H"
#include "fvMeshSubset.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class sampledInterface Declaration
\*---------------------------------------------------------------------------*/

class sampledInterface
:
    public sampledSurface
{
    // Private data

        //- Zone name/index (if restricted to zones)
        mutable cellZoneID zoneID_;

        //- For zones: patch to put exposed faces into
        mutable word exposedPatchName_;

        mutable autoPtr<interface> surfPtr_;


        // Recreated for every interface

            //- Time at last call, also track if surface needs an update
            mutable label prevTimeIndex_;

            //- Cached volfield
            mutable autoPtr<volScalarField> storedVolFieldPtr_;
            mutable const volScalarField* volFieldPtr_;

            //- Cached pointfield
            mutable const pointScalarField* pointFieldPtr_;

            // And on subsetted mesh

                //- Cached submesh
                mutable autoPtr<fvMeshSubset> subMeshPtr_;

                //- Cached volfield
                mutable autoPtr<volScalarField> storedVolSubFieldPtr_;
                mutable const volScalarField* volSubFieldPtr_;

                //- Cached pointfield
                mutable const pointScalarField* pointSubFieldPtr_;



    // Private Member Functions


        //- Create iso surface (if time has changed)
        //  Do nothing (and return false) if no update was needed
        bool updateGeometry() const;

        //- Sample volume field onto surface faces
        template<class Type>
        tmp<Field<Type>> sampleOnFaces
        (
            const interpolation<Type>& sampler
        ) const;

        //- Interpolate volume field onto surface points
        template<class Type>
        tmp<Field<Type>> sampleOnPoints
        (
            const interpolation<Type>& interpolator
        ) const;


public:

    //- Runtime type information
    TypeName("sampledInterface");


    // Constructors

        //- Construct from dictionary
        sampledInterface
        (
            const word& name,
            const polyMesh& mesh,
            const dictionary& dict
        );


    //- Destructor
    virtual ~sampledInterface();


    // Member Functions

        const interface& surface() const
        {
            return surfPtr_();
        }

        //- Does the surface need an update?
        virtual bool needsUpdate() const;

        //- Mark the surface as needing an update.
        //  May also free up unneeded data.
        //  Return false if surface was already marked as expired.
        virtual bool expire();

        //- Update the surface as required.
        //  Do nothing (and return false) if no update was needed
        virtual bool update();


        //- Points of surface
        virtual const pointField& points() const
        {
            return surface().points();
        }

        //- Faces of surface
        virtual const faceList& faces() const
        {
            return surface().surfFaces();
        }

        //- Const access to per-face zone/region information
        virtual const labelList& zoneIds() const
        {
            return Foam::emptyLabelList;
        }

        //- Face area magnitudes
        virtual const vectorField& Sf() const
        {
            return surface().Sf();
        }

        //- Face area magnitudes
        virtual const scalarField& magSf() const
        {
            return surface().magSf();
        }

        //- Face centres
        virtual const vectorField& Cf() const
        {
            return surface().Cf();
        }


        //- Lookup or read isoField.
        //  Sets volFieldPtr_ and pointFieldPtr_.
        void getIsoField();


        // Sample

            //- Sample volume field onto surface faces
            virtual tmp<scalarField> sample
            (
                const interpolation<scalar>& sampler
            ) const;

            //- Sample volume field onto surface faces
            virtual tmp<vectorField> sample
            (
                const interpolation<vector>& sampler
            ) const;

            //- Sample volume field onto surface faces
            virtual tmp<sphericalTensorField> sample
            (
                const interpolation<sphericalTensor>& sampler
            ) const;

            //- Sample volume field onto surface faces
            virtual tmp<symmTensorField> sample
            (
                const interpolation<symmTensor>& sampler
            ) const;

            //- Sample volume field onto surface faces
            virtual tmp<tensorField> sample
            (
                const interpolation<tensor>& sampler
            ) const;


        // Interpolate

            //- Interpolate volume field onto surface points
            virtual tmp<scalarField> interpolate
            (
                const interpolation<scalar>& interpolator
            ) const;

            //- Interpolate volume field onto surface points
            virtual tmp<vectorField> interpolate
            (
                const interpolation<vector>& interpolator
            ) const;

            //- Interpolate volume field onto surface points
            virtual tmp<sphericalTensorField> interpolate
            (
                const interpolation<sphericalTensor>& interpolator
            ) const;

            //- Interpolate volume field onto surface points
            virtual tmp<symmTensorField> interpolate
            (
                const interpolation<symmTensor>& interpolator
            ) const;

            //- Interpolate volume field onto surface points
            virtual tmp<tensorField> interpolate
            (
                const interpolation<tensor>& interpolator
            ) const;


    // Output

        //- Write
        virtual void print(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "sampledInterfaceTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
