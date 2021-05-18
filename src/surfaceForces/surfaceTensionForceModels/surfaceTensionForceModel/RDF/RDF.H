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
    Foam::RDF

Description
    reconstructs a signed distance function for the interface information and
    use this function to compute the curvature

SourceFiles
    RDF.C

\*---------------------------------------------------------------------------*/

#ifndef RDF_H
#define RDF_H

#include "typeInfo.H"
#include "volFields.H"
#include "dimensionedScalar.H"
#include "autoPtr.H"
#include "DynamicField.H"

#include "surfaceTensionForceModel.H"
#include "reconstructedDistanceFunction.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class RDF Declaration
\*---------------------------------------------------------------------------*/

class RDF
    : public surfaceTensionForceModel
{

private:

    //- Combine operator for Lists
    template <class type> class appendOp
    {
      public:
        void operator()(List<type>& x, const List<type>& y) const
        {
            if (y.size())
            {
                if (x.empty())
                {
                    x = y;
                }
                else
                {
                    x.append(y);
                }
            }
        }
    };

    //-
        Switch curvFromTr_;
        //- Stabilisation for normalisation of the interface normal
        const dimensionedScalar deltaN_;

        DynamicField<vector> FaceCentreField_;
        DynamicField<vector> FaceNormalField_;

        // calculates the RDF function
        reconstructedDistanceFunction& RDF_;

    // Private Member Functions
        template <typename T>
        void distributeField(const labelList& neiProcs,Field<T>& field);

        scalar distanceToSurfacePlane(const point& p);

        label closestDistToSurface(const point& p);

        labelList getNeibourProcs(const boolList& nextToInterface);

        virtual void correctContactAngle
        (
            surfaceVectorField::Boundary& nHatb,
            surfaceVectorField::Boundary& gradAlphaf
        );

        //- Disallow copy construct
        RDF(const RDF&);

        //- Disallow default bitwise assignment
        void operator=(const RDF&);

public:

    //- Runtime type information
    TypeName("RDF");

    // Constructors

        //- Construct from components
        RDF
        (
            const dictionary& dict,
            const volScalarField& alpha1,
            const surfaceScalarField& phi,
            const volVectorField& U
        );

    //- Destructor
    virtual ~RDF()
    {}

    virtual void correct();

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
