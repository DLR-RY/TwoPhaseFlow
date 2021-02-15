/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2019-2019 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
                            | Copyright (C) 2019-2019 DLR
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
    Foam::zoneDistribute

Description
    Class for parallel communication in a narrow band. It either provides a Map
    with the neighbouring values of the selected region or returns a Map of the 
    required values in global addressing. Also holds a reference to the stencil
    Before the data transfer the communation has to be set up:
    exchangeFields_.setUpCommforZone(interfaceCell_); 
    Is used in the plicRDF

    Original code supplied by Henning Scheufler, DLR (2019)

SourceFiles
    zoneDistribute.C

\*---------------------------------------------------------------------------*/

#ifndef zoneDistribute_H
#define zoneDistribute_H

#include "fvMesh.H"
#include "globalIndex.H"
#include "volFields.H"

#include "zoneCPCStencil.H"
#include "IOobject.H"
#include "MeshObject.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class zoneDistribute Declaration
\*---------------------------------------------------------------------------*/

class zoneDistribute
:   public MeshObject
    <
        fvMesh,
        TopologicalMeshObject,
        zoneDistribute
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
        // const fvMesh& mesh_;

        //- labels of the points on coupled patches
        labelList coupledBoundaryPoints_;

        //- storage of the addressing for processor-to-processor comms
        List< labelHashSet > send_;

        //- Return patch of all coupled faces.
        autoPtr<indirectPrimitivePatch> coupledFacesPatch() const;

        zoneCPCStencil& stencil_;

        const globalIndex& gblIdx_;


        //- gives patchNumber and patchFaceNumber for a given Geometric volume field
        template <typename Type>
        Type getLocalValue
        (
            const GeometricField<Type, fvPatchField, volMesh>& phi,
            const label& localIdx
        ) const;


        //- gives patchNumber and patchFaceNumber for a given Geometric volume field
        template <typename Type>
        Type faceValue
        (
            const GeometricField<Type, fvPatchField, volMesh>& phi,
            const label& localIdx
        ) const;

        //- No copy construct
        zoneDistribute(const zoneDistribute&) = delete;
  
        //- No copy assignment
        void operator=(const zoneDistribute&) = delete;

  public:
    // Constructors
        // Declare name of the class and its debug switch
        TypeName("zoneDistribute");

        //- Construct from fvMesh
        explicit zoneDistribute(const fvMesh&);

    // Selectors
        static zoneDistribute& New(const fvMesh&);

    //- Destructor
        virtual ~zoneDistribute() = default;      

        // Member functions

        //- updates the stencil with boolList the size has to match
        //- mesh.nCells()
        void setUpCommforZone(const boolList& zone,bool updateStencil=true);

        //- updates the stencil with boolList the size has to match
        //- mesh.nCells()
        void updateStencil(const boolList& zone);

        //- return stencil reference
        const labelListList& getStencil()
        {
            return stencil_;
        }

        //- return addressing reference
        const globalIndex& globalNumbering() const
        {
            return gblIdx_;
        } 

        //- gives patchNumber and patchFaceNumber for a given Geometric volume field
        template <typename Type>
        Type getValue
        (
            const GeometricField<Type, fvPatchField, volMesh>& phi,
            const Map<Type>& valuesFromOtherProc,
            const label& gblIdx
        ) const;

        //- Returns stencil and provides a Map with globalNumbering
        template <typename Type>
        Map<Field <Type> > getFields
        (
            const boolList& zone,
            const GeometricField<Type, fvPatchField, volMesh>& phi
        );

        //- Returns stencil and provides a Map with globalNumbering
        template <typename Type>
        Map< Type > getDatafromOtherProc
        (
            const boolList& zone,
            const GeometricField<Type, fvPatchField, volMesh>& phi
        );

        // virtual void updateMesh(const mapPolyMesh& mpm);

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
#include "zoneDistributeI.H"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
