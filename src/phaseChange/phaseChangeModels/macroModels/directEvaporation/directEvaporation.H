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
    Foam::directEvaporation

Description
    computes the phasechangemass from an implicit interface gradient model and
    additionally accounts for evaporation if a specified threshold
    (superheated) is exceeded


SourceFiles
    directEvaporation.C

\*---------------------------------------------------------------------------*/

#ifndef directEvaporation_H
#define directEvaporation_H

#include "typeInfo.H"
#include "macroModel.H"
#include "volFields.H"
#include "dimensionedScalar.H"
#include "autoPtr.H"

#include "fvScalarMatrix.H"
#include "fvm.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class directEvaporation Declaration
\*---------------------------------------------------------------------------*/

class directEvaporation
    : public macroModel
{
private:

    scalar superheated_;

    volScalarField boilingHeatFlow_;

    // Private Member Functions

        //- Disallow copy construct
        directEvaporation(const directEvaporation&);

        //- Disallow default bitwise assignment
        void operator=(const directEvaporation&);

public:

    //- Runtime type information
    TypeName("directEvaporation");

    // Constructors

        //- Construct from components
        directEvaporation
        (
            const phaseModel& phase1,
            const phaseModel& phase2,
            const volScalarField& p,
            singleComponentSatProp& satModel,
            const compressibleInterPhaseTransportModel& turbModel,
            const dictionary& dict
        );


    //- Destructor
    virtual ~directEvaporation()
    {}

    virtual void TSource1(fvScalarMatrix& T1Eqn);

    virtual void TSource2(fvScalarMatrix& T2Eqn);

    virtual void energySource(volScalarField& Q);

    virtual void energySource1(volScalarField& q1);

    virtual void energySource2(volScalarField& q2);

    virtual void massSource(volScalarField& rhoSource);

    virtual void alphaSource(volScalarField& rhoSource);

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
