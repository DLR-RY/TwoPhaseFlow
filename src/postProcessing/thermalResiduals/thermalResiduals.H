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
    Foam::functionObjects::thermalResiduals

Group
    grpForcesFunctionObjects

Description
    Calculates and write the heat-flux at wall patches as the
    volScalarField field 'thermalResiduals'.

    All wall patches are included by default; to restrict the calculation to
    certain patches, use the optional 'patches' entry.

Usage
    Example of function object specification:
    \verbatim
    thermalResiduals1
    {
        type        thermalResiduals;
        libs        ("libfieldFunctionObjects.so");
        ...
        patches     (".*Wall");
    }
    \endverbatim

    Where the entries comprise:
    \table
        Property | Description                | Required   | Default value
        type     | type name: thermalResiduals    | yes        |
        patches  | list of patches to process | no         | all wall patches
    \endtable

See also
    Foam::functionObject
    Foam::functionObjects::fvMeshFunctionObject
    Foam::functionObjects::writeFile
    Foam::functionObjects::timeControl

SourceFiles
    thermalResiduals.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_thermalResiduals_H
#define functionObjects_thermalResiduals_H

#include "fvMeshFunctionObject.H"
#include "writeFile.H"
#include "volFieldsFwd.H"
#include "HashSet.H"
#include "dimensionedScalar.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                       Class thermalResiduals Declaration
\*---------------------------------------------------------------------------*/

class thermalResiduals
:
    public fvMeshFunctionObject,
    public writeFile
{

protected:

    // Protected data

    word TField_;
    word alphaField_;
    word rhoField_;

    dimensionedScalar subCooled_;


    // Protected Member Functions

        //- File header information
        virtual void writeFileHeader(Ostream& os) const;



private:

    // Private member functions

        //- Disallow default bitwise copy construct
        thermalResiduals(const thermalResiduals&);

        //- Disallow default bitwise assignment
        void operator=(const thermalResiduals&);


public:

    //- Runtime type information
    TypeName("thermalResiduals");


    // Constructors

        //- Construct from Time and dictionary
        thermalResiduals
        (
            const word& name,
            const Time& runTime,
            const dictionary&
        );


    //- Destructor
    virtual ~thermalResiduals();


    // Member Functions

        //- Read the thermalResiduals data
        virtual bool read(const dictionary&);

        //- Calculate the wall heat-flux
        virtual bool execute();

        //- Write the wall heat-flux
        virtual bool write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
