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
    Foam::isoSurface

Description

    Interface reconstruction method using an isosurface of the volume faction
    field to reconstruction the interface data, i.e. interface centre and normal
    inside cells intersected by the isosurface. The default isovalue is 0.5.

    Intended for use with algebraic VoF methods like MULES, which do not
    reconstruct the interface as part of the solution process.
    
    The advantage of an isosurface is that it is continuous across cell faces.
    The disadvantages are that it does not intersect all cells containing both
    fluids (0 < alpha < 1), and that the cut cells are not cut into subsections
    in accordance with their volume fraction.

SourceFiles
    isoSurface.C

\*---------------------------------------------------------------------------*/

#ifndef isoSurface_H
#define isoSurface_H

#include "typeInfo.H"
#include "reconstructionSchemes.H"
#include "volFields.H"
#include "dimensionedScalar.H"
#include "autoPtr.H"

#include "volPointInterpolation.H"
#include "cutCellIso.H"



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace reconstruction
{

/*---------------------------------------------------------------------------*\
                Class isoSurface Declaration
\*---------------------------------------------------------------------------*/

class isoSurface
    : public reconstructionSchemes
{

private:

        //- Reference to mesh
        const fvMesh& mesh_;

        //- VOF field interpolated to mesh points
        scalarField ap_;

        cutCellIso cutCell_;

    // Switches and tolerances. Tolerances need to go into toleranceSwitches

        //- Tolerance for search of isoFace giving specified VOF value
        scalar iso_;

private:

    // Private Member Functions

        //- Disallow copy construct
        isoSurface(const isoSurface&);

        //- Disallow default bitwise assignment
        void operator=(const isoSurface&);

public:

    //- Runtime type information
    TypeName("isoSurface");

    // Constructors

        //- Construct from components
        isoSurface
        (
            volScalarField& alpha1,
            const surfaceScalarField& phi,
            const volVectorField& U,
            const dictionary& dict
        );

        //- Destructor
        virtual ~isoSurface();

        virtual void reconstruct(bool forceUpdate = true);

        virtual void mapAlphaField()
        {
            // do nothing;
        }


};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End reconstructionSchemes Foam
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
