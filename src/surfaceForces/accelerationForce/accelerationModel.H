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
    Foam::accelerationModel

Description
    abstract base class of the acceleration model

SourceFiles
    accelerationModel.C
    newAccelerationModel.C

\*---------------------------------------------------------------------------*/

#ifndef accelerationModel_H
#define accelerationModel_H

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
                Class accelerationModel Declaration
\*---------------------------------------------------------------------------*/

class accelerationModel
{
protected:

        surfaceScalarField accf_; // ghf
        volScalarField acc_; // gh

        //- Re-calculate the interface curvature
        virtual void calculateAcc();

private:

    // Private Member Functions

        //- Disallow copy construct
        accelerationModel(const accelerationModel&);

        //- Disallow default bitwise assignment
        void operator=(const accelerationModel&);

public:

    //- Runtime type information
    TypeName("accelerationModel");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            accelerationModel,
            components,
            (
                const dictionary& dict,
                const fvMesh& mesh
            ),
            (dict,mesh)
        );

    // Selectors

        //- Return a reference to the selected phaseChange model
        static autoPtr<accelerationModel> New
        (
            const dictionary& dict,
            const fvMesh& mesh
        );

    // Constructors

        //- Construct from components
        accelerationModel
        (
            const word& type,
            const dictionary& dict,
            const fvMesh& mesh
        );

    //- Destructor
    virtual ~accelerationModel()
    {}

    // Member Functions

        const surfaceScalarField& accf() const
        {
            return accf_;
        }

        const volScalarField& acc() const
        {
            return acc_;
        }

        void correct()
        {
            calculateAcc();
        }

        virtual tmp<surfaceScalarField> accelerationForce() = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
