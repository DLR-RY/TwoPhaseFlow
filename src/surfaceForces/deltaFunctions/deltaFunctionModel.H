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
    Foam::deltaFunctionModel

Description

SourceFiles
    deltaFunctionModel.C
    newDeltaFunctionModel.C

\*---------------------------------------------------------------------------*/

#ifndef deltaFunctionModel_H
#define deltaFunctionModel_H

#include "typeInfo.H"
#include "runTimeSelectionTables.H"

#include "dimensionedScalar.H"
#include "autoPtr.H"


#include "IOdictionary.H"
#include "volFields.H"
#include "surfaceFields.H"

#include "fvcSnGrad.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class deltaFunctionModel Declaration
\*---------------------------------------------------------------------------*/

class deltaFunctionModel
{
protected:

        // reference to volume fraction field
        const volScalarField& alpha1_;

        //- delta function on the faces
        surfaceScalarField deltaFunction_;

private:

    // Private Member Functions

        //- Disallow copy construct
        deltaFunctionModel(const deltaFunctionModel&);

        //- Disallow default bitwise assignment
        void operator=(const deltaFunctionModel&);

protected:

        //- Re-calculate the interface curvature
        virtual void calculateDeltaFunction();

public:

    //- Runtime type information
    TypeName("deltaFunctionModel");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            deltaFunctionModel,
            components,
            (
                const dictionary& dict,
                const volScalarField& alpha1
            ),
            (dict,alpha1)
        );


    // Selectors

        //- Return a reference to the selected phaseChange model
        static autoPtr<deltaFunctionModel> New
        (
            const dictionary& dict,
            const volScalarField& alpha1
        );


    // Constructors

        //- Construct from components
        deltaFunctionModel
        (
            const word& type,
            const dictionary& dict,
            const volScalarField& alpha1
        );


    //- Destructor
    virtual ~deltaFunctionModel()
    {}


    // Member Functions

        const surfaceScalarField& deltaFunction() const
        {
            return deltaFunction_;
        }

        tmp<surfaceScalarField> deltaFunction(const volScalarField& phi) const
        {
            return fvc::snGrad(phi);
        }

        void correct()
        {
            calculateDeltaFunction();
        }

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
