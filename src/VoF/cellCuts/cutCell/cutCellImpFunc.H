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
    Foam::cutCellImpFunc

Description


Author
    Johan Roenby, DHI, all rights reserved.

SourceFiles
    cutCellImpFunc.C

\*---------------------------------------------------------------------------*/

#ifndef cutCellImpFunc_H
#define cutCellImpFunc_H

#include "OFstream.H"
#include "cutCell.H"
#include "cutFaceImpFunc.H"
#include "fvMesh.H"
#include "scalarMatrices.H"
#include "surfaceFields.H"
#include "volFields.H"

#include "triSurface.H"
#include "triSurfaceTools.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class cutCellImpFunc Declaration
\*---------------------------------------------------------------------------*/

class cutCellImpFunc
:
    public cutCell
{
    // Private data

        //- Mesh whose cells and faces to cut at their intersection with an
        //  isosurface.
        const fvMesh& mesh_;

        //- Cell to cut
        label cellI_;

        //- Isofunction values at mesh points. f_size() = mesh_.nPoints().
        scalarField& f_;

        // implicit function
        implicitFunction& func_;

        //- Isovalue used to cut cell
        scalar cutValue_;

        //- An cutFaceIso object to get access to its face cutting functionality
        cutFaceImpFunc cutFace_;

        //- List of face centres for CutFaces
        DynamicList<point> cutFaceCentres_;

        //- List of face area vectors for isoCutFaces
        DynamicList<vector> cutFaceAreas_;

        //- Storage for subFace edges belonging to isoFace
        DynamicList<DynamicList<point>> isoFaceEdges_;

        //- Points constituting the cell-isosurface intersection (isoface)
        DynamicList<point> facePoints_;

        //- Face centre of the cutFace
        point faceCentre_;

        //- Face normal of the isoface by convention pointing from high to low
        //  values (i.e. opposite of the gradient vector).
        vector faceArea_;

        //- Cell centre of the subcell of celli which is "fully submerged", i.e.
        //  where the function value is higher than the isoValue_
        point subCellCentre_;

        //- Volume of fully submerged subcell
        scalar subCellVolume_;

        //- Volume of Fluid for cellI (subCellVolume_/mesh_.V()[cellI])
        scalar VOF_;

        //- A cell status label taking one of the values:
        //
        //  - -1: cell is fully below the isosurface
        //  -  0: cell is cut
        //  - +1: cell is fully above the isosurface
        label cellStatus_;

        triSurface freeSurf_;

    // Private Member Functions

       point bisection(point p0, point p1);

       point newtonMethod(const point& p0);

       point newtonMethod(const tensor& t,const point& p0);

       triSurface refineFreeSurfaces(triSurface& triF);

       DynamicList< point> refineCuttedFace
       (
           const vector& n,
           DynamicList< point>& cuttedFace,
           DynamicList< point>& interfacePoints
        );

        face sortInterfacePoints(const DynamicList< point>& interfacePoints,vector n);

        void getBase(const vector& n, vector& e0, vector& e1) const;

  public:
    // Constructors

        //- Construct from fvMesh and a scalarField
        //  Length of scalarField should equal number of mesh points
        cutCellImpFunc(const fvMesh&, scalarField& f,  implicitFunction& func);

        // Member functions

        //- sets internal values and returns face status
        label calcSubCell(const label& cellI, const scalar cutValue);

        //- returns subCellCentre
        point subCellCentre();

        //- returns subCellVolume
        scalar subCellVolume();

        //- returns the points of the cutting isoface
        DynamicList<point> facePoints();

        //- returns the centre of the cutting isoface
        point faceCentre();

        //- returns the area normal vector of the cutting isoface
        vector faceArea();

        //- returns cellStatus
        label cellStatus();

        //- returns volume of fluid value
        scalar VolumeOfFluid();

        //- returns cutValue
        scalar cutValue() const;

        const triSurface& triFreeSurf();

        //- resets internal values
        void clearStorage();

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
