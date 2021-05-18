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
    Foam::constantCurvature

Description
    defines a constant curvature to test if pressure velocity couplding is
    well-balanced

SourceFiles
    constantCurvature.C

\*---------------------------------------------------------------------------*/

#ifndef constantCurvature_H
#define constantCurvature_H

#include "typeInfo.H"
#include "volFields.H"
#include "dimensionedScalar.H"
#include "autoPtr.H"




#include "surfaceTensionForceModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class constantCurvature Declaration
\*---------------------------------------------------------------------------*/

class constantCurvature
    : public surfaceTensionForceModel
{

private:

    //-
        //- Stabilisation for normalisation of the interface normal
        const dimensionedScalar deltaN_;

        const scalar curv_;

    // Private Member Functions

        virtual void correctContactAngle
        (
            surfaceVectorField::Boundary& nHat,
            surfaceVectorField::Boundary& constantCurvaturef
        );


        //- Disallow copy construct
        constantCurvature(const constantCurvature&);

        //- Disallow default bitwise assignment
        void operator=(const constantCurvature&);

public:

    //- Runtime type information
    TypeName("constantCurvature");


    // Constructors

        //- Construct from components
        constantCurvature
        (
            const dictionary& dict,
            const volScalarField& alpha1,
            const surfaceScalarField& phi,
            const volVectorField& U
        );


    //- Destructor
    virtual ~constantCurvature()
    {}

    virtual void correct();

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
