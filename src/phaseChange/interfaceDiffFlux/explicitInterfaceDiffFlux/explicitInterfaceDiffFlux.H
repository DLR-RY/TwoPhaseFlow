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
    Foam::explicitInterfaceDiffFlux

Description
    computes the diffusive flux at the interface and returns a explicit source
    term

SourceFiles
    explicitInterfaceDiffFlux.C

\*---------------------------------------------------------------------------*/

#ifndef explicitInterfaceGrad_H
#define explicitInterfaceGrad_H

#include "typeInfo.H"
#include "energySourceTermModel.H"
#include "volFields.H"
#include "dimensionedScalar.H"
#include "autoPtr.H"

#include "fvScalarMatrix.H"
#include "fvm.H"

#include "zoneDistribute.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class explicitInterfaceDiffFlux Declaration
\*---------------------------------------------------------------------------*/

class explicitInterfaceDiffFlux
{

private:

    const fvMesh& mesh_;

    // Private Member Functions

        //- Disallow copy construct
        explicitInterfaceDiffFlux(const explicitInterfaceDiffFlux&);

        //- Disallow default bitwise assignment
        void operator=(const explicitInterfaceDiffFlux&);

public:

    // Constructors

        //- Construct from components
        explicitInterfaceDiffFlux
        (
            const fvMesh& mesh
        );

    //- Destructor
    virtual ~explicitInterfaceDiffFlux()
    {}

    // Member Functions
        // - return the diffusiveFlux
        tmp<volScalarField> diffusiveFlux
        (
            const boolList& markedCells,
            const volVectorField& normal,
            const volVectorField& centre,
            const volScalarField& phi,
            const volScalarField& gamma,
            const volScalarField& boundaryValue,
            const bool otherSide
        );

        // - return the diffusiveFlux
        tmp<fvScalarMatrix> diffusiveFlux
        (
            const boolList& markedCells,
            const volVectorField& normal,
            const volVectorField& centre,
            const volScalarField& phi,
            const volScalarField& gamma,
            const volScalarField& boundaryValue,
            const bool otherSide,
            const dimensionSet dim
        );

        // - return the diffusiveFlux
        void diffusiveFlux
        (
            const boolList& markedCells,
            const volVectorField& normal,
            const volVectorField& centre,
            const volScalarField& phi,
            const volScalarField& gamma,
            const volScalarField& boundaryValue,
            const bool otherSide,
            volScalarField& diffFlux,
            fvScalarMatrix& phiSource
        );

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
