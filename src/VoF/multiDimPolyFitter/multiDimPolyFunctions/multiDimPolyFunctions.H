/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2019-2019 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
                            | Copyright (C) 2019-2019 DLR
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
    Foam::multiDimPolyFunctions

Description
    base class for polynomial functions

    Original code supplied by Henning Scheufler, DLR (2019)

SourceFiles
    multiDimPolyFunctions.C

\*---------------------------------------------------------------------------*/

#ifndef multiDimPolyFunctions_H
#define multiDimPolyFunctions_H


#include "FixedList.H"
#include "vector.H"
#include "vector2D.H"
#include "scalarField.H"
#include "runTimeSelectionTables.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class multiDimPolyFunctions Declaration
\*---------------------------------------------------------------------------*/

class multiDimPolyFunctions
{

protected:
    // Private data

    // number of polynominal term
    label nTerms_;

    scalarField coeffs_;

    const Vector<label> geomDir_;

    vector geomCorrection_;

    scalarField termValues_;



public:

    //- Runtime type information
    TypeName("multiDimPolyFunctions");

    // Declare run-time constructor selection table

    // For the dictionary constructor
    declareRunTimeSelectionTable
    (
        autoPtr,
        multiDimPolyFunctions,
        word,
        (
            const Vector<label> dirs
        ),
        (dirs)
    );


 


    // Constructors

    multiDimPolyFunctions(const Vector<label> dirs);

     // Selectors
    //- Return a reference to the selected multiDimPolyFunctions
    static autoPtr<multiDimPolyFunctions> New
    (
        const word& multiDimPolyFunctionsType,
        const Vector<label> dirs
    );


    //- Destructor
    virtual ~multiDimPolyFunctions();

    virtual scalar value(const vector& vec) = 0;

    virtual const scalarField& termValues(const vector& vec) = 0;

    label& nTerms()
    {
        return nTerms_;
    }

    scalarField& coeffs()
    {
        return coeffs_;
    }

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
