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
    Foam::reconstruction::isoAlpha

Description
    Reconstructs an interface (centre and normal vectors) consisting of isosurfaces
    to match the internal fluid distribution in cells. The point values (ap_)
    are estimated by inverse distance interpolation of the VoF (alpha) field.

    Reference:
    \verbatim
        Roenby, J., Bredmose, H. and Jasak, H. (2016).
        A computational method for sharp interface advection
        Royal Society Open Science, 3
        doi 10.1098/rsos.160405
    \endverbatim

    Original code supplied by Johan Roenby, DHI (2016)

Author
    Johan Roenby, DHI, all rights reserved.
    Modified Henning Scheufler, DLR

SourceFiles
    isoAlpha.C

\*---------------------------------------------------------------------------*/

#ifndef isoAlpha_H
#define isoAlpha_H

#include "autoPtr.H"
#include "dimensionedScalar.H"
#include "reconstructionSchemes.H"
#include "typeInfo.H"
#include "volFields.H"

#include "surfaceIteratorIso.H"
#include "volPointInterpolation.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace reconstruction
{

/*---------------------------------------------------------------------------*\
                          Class isoAlpha Declaration
\*---------------------------------------------------------------------------*/

class isoAlpha
:
    public reconstructionSchemes
{
    // Private Data

        //- Reference to mesh
        const fvMesh& mesh_;

        //- VOF field interpolated to mesh points
        scalarField ap_;

        // Switches and tolerances. Tolerances need to go into toleranceSwitches

        //- Tolerance for search of isoFace giving specified VOF value
        scalar isoFaceTol_;

        //- Tolerance for marking of surface cells:
        //  Those with surfCellTol_ < alpha1 < 1 - surfCellTol_
        scalar surfCellTol_;

        //- surfaceIterator finds the isovalue for specified VOF value
        surfaceIteratorIso sIterIso_;


    // Private Member Functions

        //- No copy construct
        isoAlpha(const isoAlpha&) = delete;

        //- No copy assignment
        void operator=(const isoAlpha&) = delete;


public:

    //- Runtime type information
    TypeName("isoAlpha");

    //- Construct from components
    isoAlpha
    (
        volScalarField& alpha1,
        const surfaceScalarField& phi,
        const volVectorField& U,
        const dictionary& dict
    );


    //- Destructor
    virtual ~isoAlpha() = default;


    // Member Functions

        //- Reconstructs the interface
        virtual void reconstruct(bool forceUpdate = true);

        //- map VoF Field in case of refinement
        virtual void mapAlphaField() const
        {
            // do nothing;
        }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace reconstruction
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
