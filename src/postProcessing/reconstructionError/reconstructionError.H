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
    Foam::functionObjects::reconstructionError

Description

Author
    Henning Scheufler, DLR, all rights reserved.
    Johan Roenby, DHI, all rights reserved.

SourceFiles
    reconstructionError.C

\*---------------------------------------------------------------------------*/

#ifndef reconstructionError_H
#define reconstructionError_H


#include "writeFile.H"
#include "volFieldsFwd.H"
#include "HashSet.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class objectRegistry;

/*---------------------------------------------------------------------------*\
                       Class reconstructionError Declaration
\*---------------------------------------------------------------------------*/

class reconstructionError
:
    public functionObjects::writeFile
{
        const fvMesh& mesh_;

        scalar LNormalInf_;

        scalar LNormal1_;

        scalar LCentreInf_;

        scalar LCentre1_;

        scalar LCurvInf_;

        scalar LCurv1_;

        DynamicList<scalar> avgLNormal_;
        DynamicList<scalar> avgLCentre_;
        DynamicList<scalar> avgLCurv_;

    // Protected Member Functions

        //- File header information
        void writeFileHeader(Ostream& os) const;



private:

    // Private member functions

        //- Disallow default bitwise copy construct
        reconstructionError(const reconstructionError&);

        //- Disallow default bitwise assignment
        void operator=(const reconstructionError&);



public:


    // Constructors

        //- Construct from Time and dictionary
        reconstructionError
        (
            const fvMesh& mesh,
            const objectRegistry& obr,
            const dictionary& dict
        );


    //- Destructor
    virtual ~reconstructionError();


    // Member Functions


        //- Calculate the wall heat-flux
        void calcError
        (
            const dictionary& dict,
            const volVectorField& centre,
            const volVectorField& normal,
            const volScalarField& curv,
            bool randomMode
        );

        //- Write the wall heat-flux
        void write();

        //- Write the wall heat-flux
        void write
        (
            const scalar k1,
            const scalar k2
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
