/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2017 DHI
    Copyright (C) 2018-2019 Johan Roenby
    Copyright (C) 2019-2020 DLR
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
    Foam::cutFaceIso

Description
    Class for cutting a face, faceI, of an fvMesh, mesh_, at its intersection
    with an isosurface defined by the mesh point values f_ and the isovalue,
    cutValue.

    Reference:
    \verbatim
        Roenby, J., Bredmose, H. and Jasak, H. (2016).
        A computational method for sharp interface advection
        Royal Society Open Science, 3
        doi 10.1098/rsos.160405

    \endverbatim

    Original code supplied by Johan Roenby, DHI (2016)

SourceFiles
    cutFaceIso.C

\*---------------------------------------------------------------------------*/

#ifndef cutFaceIso_H
#define cutFaceIso_H

#include "fvMesh.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "cutFace.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class cutFaceIso Declaration
\*---------------------------------------------------------------------------*/

class cutFaceIso
:
    public cutFace
{
    // Private Data

        //- Mesh whose cells and faces to cut at their intersection with an
        //- isoface
        const fvMesh& mesh_;

        //- Isofunction values at mesh points. f_size() = mesh_.nPoints().
        scalarField& f_;

        //- Storage for centre of subface
        point subFaceCentre_;

        //- Storage for area vector of subface
        vector subFaceArea_;

        //- Storage for subFacePoints
        DynamicList<point> subFacePoints_;

        //- Storage for subFacePoints
        DynamicList<point> surfacePoints_;

        //- Storage for pointStatus_ cuts the cell at 0
        DynamicList<scalar> pointStatus_;

        //- Storage of the edge weight
        DynamicList<scalar> weight_;

        //- A face status label taking one of the values:
        //
        //   -1: face is fully below the isosurface
        //    0: face is cut, i.e. has values larger and smaller than isoValue_
        //   +1: face is fully above the isosurface
        label faceStatus_;


public:

    // Constructors

        //- Construct from fvMesh and a scalarField
        //  length of scalarField should equal number of mesh points
        cutFaceIso(const fvMesh& mesh, scalarField& f);


    // Member Functions

        //- Calculate cut points along edges of faceI
        label calcSubFace
        (
            const label faceI,
            const scalar cutValue
        );

        //- Returns centre of cutted face
        const point& subFaceCentre() const;

        //- Returns area vector of cutted face
        const vector& subFaceArea() const;

        //- Returns the cut edge of the cutted face
        const DynamicList<point>& subFacePoints() const;

        //- Returns point of the face in sorted of cutted face
        const DynamicList<point>& surfacePoints() const;

        //- Resets internal variables
        void clearStorage();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
