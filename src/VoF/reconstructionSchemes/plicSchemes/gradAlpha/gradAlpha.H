/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
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
    Foam::reconstruction::gradAlpha

Description
    Reconstructs an interface (centre and normal vector) consisting of planes
    to match the internal fluid distribution in cells. The interface normals
    are estimated by least square gradient scheme on the VoF Field (alpha).
    Also known as Young method

    Original code supplied by Henning Scheufler, DLR (2019)

SourceFiles
    gradAlpha.C

\*---------------------------------------------------------------------------*/

#ifndef gradAlpha_H
#define gradAlpha_H

#include "typeInfo.H"
#include "reconstructionSchemes.H"
#include "volFields.H"
#include "dimensionedScalar.H"
#include "autoPtr.H"
#include "surfaceIteratorPLIC.H"
#include "zoneDistribute.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace reconstruction
{

/*---------------------------------------------------------------------------*\
                          Class gradAlpha Declaration
\*---------------------------------------------------------------------------*/

class gradAlpha
:
    public reconstructionSchemes
{
    // Private Data

        //- Reference to mesh
        const fvMesh& mesh_;

        //- Interpolation object from cell centres to points
        DynamicField<vector> interfaceNormal_;


    // Switches and tolerances. Tolerances need to go into toleranceSwitches

        //- Tolerance for search of isoFace giving specified VOF value
        scalar isoFaceTol_;

        //- Tolerance for marking of surface cells:
        //  Those with surfCellTol_ < alpha1 < 1 - surfCellTol_
        scalar surfCellTol_;

        //- SurfaceIterator finds the plane centre for specified VOF value
        surfaceIteratorPLIC sIterPLIC_;

        //- Compute gradient at the surfaces
        void gradSurf(const volScalarField& phi);


    // Private Member Functions

        //- No copy construct
        gradAlpha(const gradAlpha&) = delete;

        //- No copy assignment
        void operator=(const gradAlpha&) = delete;


public:

    //- Runtime type information
    TypeName("gradAlpha");

    //- Construct from components
    gradAlpha
    (
        volScalarField& alpha1,
        const surfaceScalarField& phi,
        const volVectorField& U,
        const dictionary& dict
    );

    //- Destructor
    virtual ~gradAlpha() = default;


    // Member Functions

        //- Reconstruct interface
        virtual void reconstruct(bool forceUpdate = true);

        //- map VoF Field in case of refinement
        virtual void mapAlphaField() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace reconstruction
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
