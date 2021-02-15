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
    Foam::zoneDistributePoints

Description
    Class for parallel communication in a narrow band for pointValues

Author
    Henning Scheufler, DLR, all rights reserved.

SourceFiles
    zoneDistributePoints.C

\*---------------------------------------------------------------------------*/

#ifndef zoneDistributePoints_H
#define zoneDistributePoints_H

#include "fvMesh.H"
#include "volFields.H"

#include "IOobject.H"
#include "MeshObject.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class zoneDistributePoints Declaration
\*---------------------------------------------------------------------------*/

class zoneDistributePoints
:   public MeshObject
    <
        fvMesh,
        TopologicalMeshObject,
        zoneDistributePoints
    >
{

  private:

    // Private data

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

        //- Reference to mesh
        const fvMesh& mesh_;

        //- boundaryPoints without empty and wedge points
        labelList boundaryPoints_;

        //- is point and patch empty and wege are excluded
        boolList isBoundaryPoint_;

        //- valid boundary face excludes processor wedge and empty
        //  size is nFaces - internalFaces
        boolList validBoundaryFace_;

        //- Return patch with non-empty and non-wedge faces.
        autoPtr<indirectPrimitivePatch> nonEmptyWedgePatch() const;

        //- return faceValue
        template <typename Type>
        Type faceValue
        (
            const GeometricField<Type, fvPatchField, volMesh>& phi,
            const label faceI
        ) const;



  public:
    // Constructors
        //- Declare name of the class and its debug switch
        TypeName("zoneDistributePoints");

        //- Construct from fvMesh
        zoneDistributePoints(const fvMesh&);

        //- return boundaryPoints without wedge and empty patchPoints
        const boolList& isBoundaryPoint()
        {
            return isBoundaryPoint_;
        }


        //- Returns Map with cellPointValues
        template <typename Type>
        Map<Field <Type> > getPointDatafromOtherProc
        (
            const boolList& zone,
            const GeometricField<Type, fvPatchField, volMesh>& phi
        );

        // //- update mesh in case of topochange
        // virtual void updateMesh(const mapPolyMesh& mpm);

        // //- movePoints
        // virtual bool movePoints()
        // {
        //     // do nothing;
        //     return false;
        // }






};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#include "zoneDistributePointsI.H"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
