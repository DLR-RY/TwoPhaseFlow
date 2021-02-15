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
    Foam::implicitGrad

Description
    computes the phasechangemass from an implicit interface gradient model

SourceFiles
    implicitGrad.C

\*---------------------------------------------------------------------------*/

#ifndef implicitGrad_H
#define implicitGrad_H

#include "typeInfo.H"
#include "energySourceTermModel.H"
#include "volFields.H"
#include "dimensionedScalar.H"
#include "autoPtr.H"

#include "fvScalarMatrix.H"
#include "fvm.H"

#include "reconstructionSchemes.H"
#include "zoneDistribute.H"
#include "implicitInterfaceDiffFlux.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class implicitGrad Declaration
\*---------------------------------------------------------------------------*/

class implicitGrad
    : public energySourceTermModel
{

private:

    implicitInterfaceDiffFlux impDiffFlux_;

    // Private Member Functions

        //- Disallow copy construct
        implicitGrad(const implicitGrad&);

        //- Disallow default bitwise assignment
        void operator=(const implicitGrad&);


public:

    //- Runtime type information
    TypeName("implicitGrad");

    // Constructors

        //- Construct from components
        implicitGrad
        (
            const phaseModel& phase1,
            const phaseModel& phase2,
            const compressibleInterPhaseTransportModel& turbModel,
            const volScalarField& p,
            singleComponentSatProp& satModel,
            reconstructionSchemes& surf,
            const dictionary& dict
        );


    //- Destructor
    virtual ~implicitGrad()
    {}

    // Member Functions
        virtual tmp<fvScalarMatrix> TSource1();

        virtual tmp<fvScalarMatrix> TSource2();

        virtual tmp<volScalarField> energySource();

        virtual tmp<volScalarField> energyFlux1();

        virtual tmp<volScalarField> energyFlux2();

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
