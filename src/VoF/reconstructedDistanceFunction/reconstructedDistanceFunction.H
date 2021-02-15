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
    Foam::reconstructedDistanceFunction

Description
    Calculates a reconstructed distance function if the cells a more than
    two cell layer a ways the octree based method needs to be used. Also
    marks the cell layers near the interface cells.

Author
    Henning Scheufler, DLR, all rights reserved.

SourceFiles
    reconstructedDistanceFunction.C

\*---------------------------------------------------------------------------*/

#ifndef reconstructedDistanceFunction_H
#define reconstructedDistanceFunction_H

#include "fvMesh.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "globalIndex.H"
#include "Map.H"
#include "zoneDistribute.H"
#include "dimensionedScalar.H"
#include "typeInfo.H"



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class reconstructedDistanceFunction Declaration
\*---------------------------------------------------------------------------*/

class reconstructedDistanceFunction
: public volScalarField
{


private:

    // Private data

    // reference to mesh
    const fvMesh& mesh_;

    // stores the coupled boundary points which have to be synced
    labelList coupledBoundaryPoints_;

    volScalarField cellDistLevel_;

    boolList nextToInterface_;

    //- Return patch of all coupled faces.
    autoPtr<indirectPrimitivePatch> coupledFacesPatch() const;





public:

    // Constructors
    // Declare name of the class and its debug switch
    ClassName("reconstructedDistanceFunction");

    //- Construct from fvMesh and a scalarField
    //  Length of scalarField should equal number of mesh points
    reconstructedDistanceFunction(const fvMesh& mesh);

    // Selectors
    static reconstructedDistanceFunction& New(const fvMesh&);


    // Member functions

    //void markCellsNearSurf(const boolList& interfaceCells,const label neiRingLevel);

    void markCellsNearSurf(const boolList& interfaceCells,const label neiRingLevel);

    const volScalarField& constructRDF
    (
        const boolList& interfaceCells,
        const volVectorField& centre,
        const volVectorField& normal,
        const label neiRingLevel,
        zoneDistribute& distribute
    );

    const volScalarField& constructRDF
    (
        const boolList& nextToInterface,
        const volVectorField& centre,
        const volVectorField& normal,
        zoneDistribute& distribute,
        bool updateStencil=true
    );


    const volScalarField& constructRDFOctree
    (
        const boolList& nextToInterface,
        const pointField& centre,
        const vectorField& normals
    );

    void updateContactAngle
    (
        const volScalarField& alpha,
        const volVectorField& U,
        surfaceVectorField::Boundary& nHatb
    );

    const volScalarField& cellDistLevel()
    {
        return cellDistLevel_;
    }

    const boolList& nextToInterface()
    {
        return nextToInterface_;
    }

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
