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
    Foam::curvatureModel

Description

SourceFiles
    curvatureModel.C
    newCurvatureModel.C

\*---------------------------------------------------------------------------*/

#ifndef curvatureModel_H
#define curvatureModel_H

#include "typeInfo.H"
#include "runTimeSelectionTables.H"
#include "volFields.H"
#include "dimensionedScalar.H"
#include "autoPtr.H"
#include "Pair.H"

#include "IOdictionary.H"
#include "volFields.H"
#include "surfaceFields.H"



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class curvatureModel Declaration
\*---------------------------------------------------------------------------*/

class curvatureModel
{

protected:

        //- reference to volume fraction field
        const volScalarField& alpha1_;

        //- reference to face velocity
        const surfaceScalarField& phi_;

        //- reference to velcoity
        const volVectorField& U_;

        //- surface face normal on the face
        surfaceScalarField nHatf_;

        //- curvature field
        volScalarField K_;

        //- curvature field on faces
        surfaceScalarField Kf_;

private:

    // Private Member Functions

        //- Disallow copy construct
        curvatureModel(const curvatureModel&);

        //- Disallow default bitwise assignment
        void operator=(const curvatureModel&);

protected:

        //- Re-calculate the interface curvature
        virtual void calculateK();

public:

    //- Runtime type information
    TypeName("curvatureModel");

    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            curvatureModel,
            components,
            (
                const dictionary& dict,
                const volScalarField& alpha1,
                const surfaceScalarField& phi,
                const volVectorField& U
            ),
            (dict,alpha1,phi,U)
        );

    // Selectors

        //- Return a reference to the selected phaseChange model
        static autoPtr<curvatureModel> New
        (
            const dictionary& dict,
            const volScalarField& alpha1,
            const surfaceScalarField& phi,
            const volVectorField& U
        );

    // Constructors

        //- Construct from components
        curvatureModel
        (
            const word& type,
            const dictionary& dict,
            const volScalarField& alpha1,
            const surfaceScalarField& phi,
            const volVectorField& U
        );

    //- Destructor
    virtual ~curvatureModel()
    {}

    // Member Functions

        const volScalarField& K() const
        {
            return K_;
        }

        const surfaceScalarField& Kf() const
        {
            return Kf_;
        }

        const surfaceScalarField& nHatf() const
        {
            return nHatf_;
        }

        void correct()
        {
            calculateK();
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
