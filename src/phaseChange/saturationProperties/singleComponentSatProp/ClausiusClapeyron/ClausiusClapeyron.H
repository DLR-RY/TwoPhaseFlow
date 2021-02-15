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
    Foam::ClausiusClapeyron

Description

    Clausius Clapeyron Relation:

    ln(p) - ln(p0) = L/R *(1/T0 - 1/T) // T0 and p0 are constants
    ln(p) = L/(R-T0) + ln(p0) + L/(R*T)
    -> ln(p) = C0 - L/(R*T) with C0 = L/(R*T0) + ln(p0)

SourceFiles
    ClausiusClapeyron.C

\*---------------------------------------------------------------------------*/

#ifndef ClausiusClapeyron_H
#define ClausiusClapeyron_H

#include "typeInfo.H"
#include "singleComponentSatProp.H"
#include "volFields.H"
#include "dimensionedScalar.H"
#include "autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class ClausiusClapeyron Declaration
\*---------------------------------------------------------------------------*/

class ClausiusClapeyron
    : public singleComponentSatProp
{

private:

    //-
    scalar pSat0_;

    //
    scalar TSat0_;

    //-
    scalar L0_;

    //-
    scalar R_;

    scalar C0_;

    scalar C0();




    // Private Member Functions

        //- Disallow copy construct
        ClausiusClapeyron(const ClausiusClapeyron&);

        //- Disallow default bitwise assignment
        void operator=(const ClausiusClapeyron&);


public:

    //- Runtime type information
    TypeName("ClausiusClapeyron");


    // Constructors

        //- Construct from components
        ClausiusClapeyron
        (
            const fvMesh& mesh,
            const dictionary& dict
        );


    //- Destructor
    virtual ~ClausiusClapeyron()
    {}


    // Member Functions

        //- Read the transportProperties dictionary and update
        //-
        virtual scalar satP (scalar T) const;

        //-
        virtual scalar satT (scalar p) const;

        //
        virtual scalar satL (scalar p) const;


        // virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
