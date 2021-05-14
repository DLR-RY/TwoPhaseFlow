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
    Foam::reconstructionSchemes

Description
    Original code supplied by Henning Scheufler, DLR (2019)

SourceFiles
    reconstructionSchemes.C

\*---------------------------------------------------------------------------*/

#ifndef reconstructionSchemes_H
#define reconstructionSchemes_H

#include "typeInfo.H"
#include "runTimeSelectionTables.H"
#include "dimensionedScalar.H"
#include "volFields.H"
#include "IOdictionary.H"
#include "DynamicField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class reconstructionSchemes Declaration
\*---------------------------------------------------------------------------*/

class reconstructionSchemes
:
    public IOdictionary
{
    // Private Data

        //- Coefficients dictionary
        dictionary reconstructionSchemesCoeffs_;


protected:

    // Protected Data

        //- Reference to the VoF Field
        volScalarField& alpha1_;

        //- Reference to the face fluxes
        const surfaceScalarField& phi_;

        //- Reference to the velocity field
        const volVectorField& U_;

        //- Interface area normals
        volVectorField normal_;

        //- Interface centres
        volVectorField centre_;

        //- Is interface cell?
        boolList interfaceCell_;

        //- List of cell labels that have an interface
        DynamicField<label> interfaceLabels_;

        //- Store timeindex/iteration to avoid multiple reconstruction
        mutable Pair<label> timeIndexAndIter_;


    // Protected Member Functions

        //- Is the interface already up-to-date?
        bool alreadyReconstructed(bool forceUpdate = true) const;

        //- No copy construct
        reconstructionSchemes(const reconstructionSchemes&) = delete;

        //- No copy assignment
        void operator=(const reconstructionSchemes&) = delete;




public:

    //- Runtime type information
    TypeName("reconstructionSchemes");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            reconstructionSchemes,
            components,
            (
                volScalarField& alpha1,
                const surfaceScalarField& phi,
                const volVectorField& U,
                const dictionary& dict
            ),
            (alpha1, phi, U, dict)
        );


    // Selectors

        //- Return a reference to the selected phaseChange model
        static autoPtr<reconstructionSchemes> New
        (
            volScalarField& alpha1,
            const surfaceScalarField& phi,
            const volVectorField& U,
            const dictionary& dict
        );


    // Constructors

        //- Construct from components
        reconstructionSchemes
        (
            const word& type,
            volScalarField& alpha1,
            const surfaceScalarField& phi,
            const volVectorField& U,
            const dictionary& dict
        );


    //- Destructor
    virtual ~reconstructionSchemes() = default;


    // Member Functions

        //- Access to the model dictionary
        dictionary& modelDict()
        {
            return reconstructionSchemesCoeffs_;
        }

        //- Const access to the model dictionary
        const dictionary& modelDict() const
        {
            return reconstructionSchemesCoeffs_;
        }

        //- Reconstruct the interface
        virtual void reconstruct(bool forceUpdate = true) = 0;

        //- const-Reference to interface normals
        const volVectorField& normal() const
        {
            return normal_;
        }

        //- const-Reference to interface centres
        const volVectorField& centre() const
        {
            return centre_;
        }

        //- Reference to interface normals
        volVectorField& normal()
        {
            return normal_;
        }

        //- Reference to interface centres
        volVectorField& centre()
        {
            return centre_;
        }

        //- Does the cell contain interface
        const boolList& interfaceCell() const
        {
            return interfaceCell_;
        }

        //- List of cells with an interface
        const DynamicField<label>& interfaceLabels() const
        {
            return interfaceLabels_;
        }

        //- Map VoF Field in case of refinement
        virtual void mapAlphaField() const
        {}


};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
