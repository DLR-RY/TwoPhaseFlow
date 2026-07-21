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

    Water saturation curve by Roche-Magnus

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
