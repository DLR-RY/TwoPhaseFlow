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
    Foam::fitParaboloid

Description
    estimates the curvature by fitting a paraboloid in the interface centres

SourceFiles
    fitParaboloid.C

\*---------------------------------------------------------------------------*/

#ifndef fitParaboloid_H
#define fitParaboloid_H

#include "typeInfo.H"
#include "volFields.H"
#include "dimensionedScalar.H"
#include "autoPtr.H"

#include "DynamicField.H"

#include "surfaceTensionForceModel.H"
#include "zoneDistribute.H"
#include "cartesianCS.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class fitParaboloid Declaration
\*---------------------------------------------------------------------------*/

class fitParaboloid
    : public surfaceTensionForceModel
{

    private:

        //-
        //- Stabilisation for normalisation of the interface normal
        const dimensionedScalar deltaN_;

        //- update contact angle
        virtual void correctContactAngle
        (
            surfaceVectorField::Boundary& nHat,
            surfaceVectorField::Boundary& gradAlphaf
        );

        //- update contact angle
        void correctContactAngle
        (
            volVectorField& normal,
            volVectorField& centre
        );

        //- compute curvature from the fit data
        scalar calcCurvature
        (
            const scalarField& fit
        );

        //- get the faceCentre
        vectorField getFaceCentres
        (
            const globalIndex& globalNumbering,
            const DynamicList < label >& stencil,
            const volVectorField& faceCentres,
            const Map < vector >& map
        );

        //- Disallow copy construct
        fitParaboloid(const fitParaboloid&);

        //- Disallow default bitwise assignment
        void operator=(const fitParaboloid&);

    public:

        //- Runtime type information
        TypeName("fitParaboloid");

    // Constructors
        //- Construct from components
        fitParaboloid
        (
            const dictionary& dict,
            const volScalarField& alpha1,
            const surfaceScalarField& phi,
            const volVectorField& U
        );

        //- Destructor
        virtual ~fitParaboloid()
        {}


        virtual void correct();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
