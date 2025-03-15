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
    Foam::improvedGravity

Description
    constant gravitional forces

SourceFiles
    improvedGravity.C

\*---------------------------------------------------------------------------*/

#ifndef improvedGravity_H
#define improvedGravity_H

#include "typeInfo.H"
#include "volFields.H"
#include "dimensionedScalar.H"
#include "autoPtr.H"
#include "uniformDimensionedFields.H"

#include "accelerationForceModel.H"
#include "alphaCSF.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class improvedGravity Declaration
\*---------------------------------------------------------------------------*/

class improvedGravity
    : public accelerationForceModel
{

        //- const access dict
        const dictionary& improvedGravityDict_;

        //- improvedGravity vector
        dimensionedVector g_;

        //- Href Field
        uniformDimensionedScalarField hRef_;

    // Private Member Functions

        //- Disallow copy construct
        improvedGravity(const improvedGravity&);

        //- Disallow default bitwise assignment
        void operator=(const improvedGravity&);


        protected:

        //- Re-calculate the interface curvature
        virtual void calculateAcc();


public:

    //- Runtime type information
    TypeName("improvedGravity");


    // Constructors

        //- Construct from components
        improvedGravity
        (
            const dictionary& dict,
            const fvMesh& mesh
        );


    //- Destructor
    virtual ~improvedGravity()
    {}

    virtual tmp<surfaceScalarField> accelerationForce();


};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
