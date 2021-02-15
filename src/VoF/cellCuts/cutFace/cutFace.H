/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016-2017 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
                isoAdvector | Copyright (C) 2016-2017 DHI
              Modified work | Copyright (C) 2018-2019 Johan Roenby
              Modified work | Copyright (C) 2019 DLR
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
    Foam::cutFace

Description
    Base class for cutting a face, faceI, of an fvMesh, mesh_, at its
    intersections

    Original code supplied by
        Johan Roenby, DHI (2016)
        Henning Scheufler, DLR (2019)



SourceFiles
    cutFace.C

\*---------------------------------------------------------------------------*/

#ifndef cutFace_H
#define cutFace_H

#include "OFstream.H"
#include "fvMesh.H"
#include "surfaceFields.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class cutFaces Declaration
\*---------------------------------------------------------------------------*/

class cutFace
{
    // Private data

        //- Reference to mesh
        const fvMesh& mesh_;


protected:

    //- Calculate cut points along edges of face with pointStatus, pointfield
    //- and computes geometric information and the face status where:
    // -  -1: faceStatus is fully below the isosurface
    // -   0: faceStatus is cut, i.e. has values larger and smaller than isoValue
    // -  +1: faceStatus is fully above the isosurface
    void calcSubFace
    (
        const label& faceI,
        const scalarList& pointStatus,
        label firstFullySubmergedPoint,
        DynamicList<point>& subFacePoints,
        DynamicList<point>& surfacePoints,
        label& faceStatus,
        vector& subFaceCentre,
        vector& subFaceArea
     );

    //- Calculate cut points along edges of face with pointStatus, pointfield, 
    //- weights and computes geometric information and the face status where:
    // -  -1: faceStatus is fully below the isosurface
    // -   0: faceStatus is cut, i.e. has values larger and smaller than isoValue
    // -  +1: faceStatus is fully above the isosurface
    void calcSubFace
    (
        const label& faceI,
        const scalarList& pointStatus,
        const scalarList& weights,
        label firstFullySubmergedPoint,
        DynamicList<point>& subFacePoints,
        DynamicList<point>& surfacePoints,
        label& faceStatus,
        vector& subFaceCentre,
        vector& subFaceArea
    );

    //- Calculates centre and normal of the face
    void calcSubFaceCentreAndArea
    (
        DynamicList<point>& subFacePoints,
        vector& subFaceCentre,
        vector& subFaceArea
    );

    //- Calculate cut points along edges of face with pointStatus, pointfield
    //- and computes geometric information and the face status where:
    // -  -1: faceStatus is fully below the isosurface
    // -   0: faceStatus is cut, i.e. has values larger and smaller than isoValue
    // -  +1: faceStatus is fully above the isosurface
    void calcSubFace
    (
        const face& f,
        const pointField& points,
        const scalarList& pointStatus,
        label firstFullySubmergedPoint,
        DynamicList<point>& subFacePoints,
        DynamicList<point>& surfacePoints,
        label& faceStatus,
        vector& subFaceCentre,
        vector& subFaceArea
    );


public:

    //- Construct from fvMesh
    cutFace(const fvMesh& mesh);

    // Static data
  
    static int debug;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
