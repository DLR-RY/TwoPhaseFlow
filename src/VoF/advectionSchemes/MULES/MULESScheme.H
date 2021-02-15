/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
     \\/     M anipulation  |
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
    Foam::MULES

Description
    interface compression volume of fluid method

SourceFiles
    MULESScheme.C

\*---------------------------------------------------------------------------*/

#ifndef MULESScheme_H
#define MULESScheme_H

#include "typeInfo.H"
#include "advectionSchemes.H"
#include "volFields.H"
#include "dimensionedScalar.H"
#include "autoPtr.H"

#include "volPointInterpolation.H"
#include "reconstructionSchemes.H"



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace advection
{

/*---------------------------------------------------------------------------*\
                Class MULESScheme Declaration
\*---------------------------------------------------------------------------*/

class MULESScheme
    : public advectionSchemes
{

private:

        //- Reference to mesh
        const fvMesh& mesh_;

        volScalarField alpha2_;

        surfaceScalarField talphaPhiCorr0_;

        const dimensionedScalar deltaN_;

        surfaceScalarField nHatf_;


        //- Reference to mesh
        bool LTS_;

        label nAlphaCorr_;

        label nAlphaSubCycles_;

        bool MULESCorr_;

        // Apply the compression correction from the previous iteration
        // Improves efficiency for steady-simulations but can only be applied
        // once the alpha field is reasonably steady, i.e. fully developed
        bool alphaApplyPrevCorr_;//

        // Isotropic compression coefficient
        scalar icAlpha_;

        //- Compression coefficient
        scalar cAlpha_;

        void alphaEqn();

        void updateNHatf();


private:

    // Private Member Functions

        //- Disallow copy construct
        MULESScheme(const MULESScheme&);

        //- Disallow default bitwise assignment
        void operator=(const MULESScheme&);



public:

    //- Runtime type information
    TypeName("MULESScheme");


    // Constructors

        //- Construct from components
        MULESScheme
        (
                volScalarField &alpha1,
                const surfaceScalarField& phi,
                const volVectorField& U
        );


    //- Destructor
    virtual ~MULESScheme();



    virtual void advect();

    virtual void advect(const volScalarField::Internal& Sp,const volScalarField::Internal& Su);


};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End advection Foam
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
