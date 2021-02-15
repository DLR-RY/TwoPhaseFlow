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
    Foam::implicitInterfaceDiffFlux

Description
    computes the diffusive flux at the interface and returns a implicit source
    term

SourceFiles
    implicitInterfaceDiffFlux.C

\*---------------------------------------------------------------------------*/

#ifndef implicitInterfaceDiffFlux_H
#define implicitInterfaceDiffFlux_H

#include "volFields.H"
#include "dimensionedScalar.H"
#include "autoPtr.H"

#include "fvScalarMatrix.H"
#include "fvm.H"

#include "reconstructionSchemes.H"
#include "zoneDistribute.H"
#include "tmp.H"
#include "volFields.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class implicitInterfaceDiffFlux Declaration
\*---------------------------------------------------------------------------*/

class implicitInterfaceDiffFlux
{

private:

    const fvMesh& mesh_;

    label getCellLabel
    (
        const label& localIndex
    );

    scalar diffusiveFlux
    (
        label cellI,
        const volVectorField& centre,
        const volVectorField& normal,
        const volScalarField& phi,
        const volScalarField& gamma,
        const bool otherSide
    );

    void getNeigbours
    (
        const label celli,
        const vector faceCentre,
        const vector faceNormal,
        const volVectorField& centre,
        const volScalarField& phi,
        const Map<vector>& mapCC,
        const Map<scalar>& mapPhi,
        zoneDistribute& exchangeFields,
        DynamicField<scalar>& phiField,
        DynamicField<vector>& distField,
        DynamicField<label>& indicies
    );

    // Private Member Functions

        //- Disallow copy construct
        implicitInterfaceDiffFlux(const implicitInterfaceDiffFlux&);

        //- Disallow default bitwise assignment
        void operator=(const implicitInterfaceDiffFlux&);

public:

    // Constructors

        //- Construct from components
        implicitInterfaceDiffFlux
        (
            const fvMesh& mesh
        );

    //- Destructor
    virtual ~implicitInterfaceDiffFlux()
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
