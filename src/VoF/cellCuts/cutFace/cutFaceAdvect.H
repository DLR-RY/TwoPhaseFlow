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
    Foam::cutFaceAdvect

Description
    Calculates the face fluxes

    Reference:
    \verbatim
        Roenby, J., Bredmose, H. and Jasak, H. (2016).
        A computational method for sharp interface advection
        Royal Society Open Science, 3
        doi 10.1098/rsos.160405

    \endverbatim

    Original code supplied by Johan Roenby, DHI (2016)

SourceFiles
    cutFaceAdvect.C

\*---------------------------------------------------------------------------*/

#ifndef cutFaceAdvect_H
#define cutFaceAdvect_H

#include "cutFace.H"
#include "fvMesh.H"
#include "surfaceFields.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class cutFaceAdvect Declaration
\*---------------------------------------------------------------------------*/

class cutFaceAdvect
:
    public cutFace
{
    // Private Data

        typedef DynamicList<label> DynamicLabelList;
        typedef DynamicList<scalar> DynamicScalarList;
        typedef DynamicList<vector> DynamicVectorList;
        typedef DynamicList<point> DynamicPointList;

        //- Mesh whose cells and faces to cut at their intersection
        //- with an isoface
        const fvMesh& mesh_;

        //- Reference to VoF field
        const volScalarField& alpha1_;

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

        //- Storage for arrival Times
        DynamicList<scalar> pTimes_;

        //- A face status label taking one of the values:
        //
        //   -1: face is fully below the isosurface
        //    0: face is cut, i.e. has values larger and smaller than isoValue_
        //   +1: face is fully above the isosurface
        label faceStatus_;


    // Private Member Functions

        //- Write faces to vtk files
        void isoFacesToFile
        (
            const DynamicList<List<point>>& isoFacePts,
            const word& fileName,
            const word& fileDir
        ) const;


  public:

    // Constructors

        //- Construct from fvMesh and a scalarField
        //  Length of scalarField should equal number of mesh points
        cutFaceAdvect(const fvMesh& mesh, const volScalarField& alpha1);


    // Member Functions

        //- Calculates cut centre and cut area (plicReconstruction)
        //  \return face status
        label calcSubFace
        (
            const label faceI,
            const vector& normal,
            const vector& base
        );

        //- Calculates cut centre and cut area (iso reconstruction)
        //  \return face status
        label calcSubFace
        (
            const face& f,
            const pointField& points,
            const scalarField& val,
            const scalar cutValue
        );

        //- Calculates cut centre and cut area (iso reconstruction)
        //  \return face status
        label calcSubFace
        (
            const label faceI,
            const label sign,
            const scalar cutValue
        );

        //- Calculate time integrated flux for a face
        scalar timeIntegratedFaceFlux
        (
            const label faceI,
            const vector& x0,
            const vector& n0,
            const scalar Un0,
            const scalar dt,
            const scalar phi,
            const scalar magSf
        );

        //- Calculate time integrated flux for a face
        scalar timeIntegratedFaceFlux
        (
            const label faceI,
            const scalarField& times,
            const scalar Un0,
            const scalar dt,
            const scalar phi,
            const scalar magSf
        );

        //- Calculate time integrated area for a face
        scalar timeIntegratedArea
        (
            const label faceI,
            const scalar dt,
            const scalar magSf,
            const scalar Un0
        );

        //- Calculate time integrated area for a face
        scalar timeIntegratedArea
        (
            const pointField& fPts,
            const scalarField& pTimes,
            const scalar dt,
            const scalar magSf,
            const scalar Un0
        );

        //- For face with vertices p calculate its area and integrated area
        //  between isocutting lines with isovalues f0 and f1 given vertex
        //  values f
        void quadAreaCoeffs
        (
            const DynamicPointList& pf0,
            const DynamicPointList& pf1,
            scalar& quadArea,
            scalar& intQuadArea
        ) const;

        //- Get cutPoints of face
        void cutPoints
        (
            const label faceI,
            const scalar f0,
            DynamicList<point>& cutPoints
        );

        //- Get cutPoints of face
        void cutPoints
        (
            const pointField& pts,
            const scalarField& f,
            const scalar f0,
            DynamicList<point>& cutPoints
        );

        //- Returns centre of cutted face
        const point& subFaceCentre() const;

        //- Returns area vector of cutted face
        const vector& subFaceArea() const;

        //- Returns the cut edge of the cutted face
        const DynamicList<point>& subFacePoints() const;

        //- Returns point of the cutted face in sorted order
        const DynamicList<point>& surfacePoints() const;

        //- Resets internal variables
        void clearStorage();
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
