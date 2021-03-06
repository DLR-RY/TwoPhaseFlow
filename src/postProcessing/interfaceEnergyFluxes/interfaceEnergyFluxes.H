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
    Foam::functionObjects::interfaceEnergyFluxes

Group
    grpForcesFunctionObjects

Description
    Calculates and write the heat-flux at wall patches as the
    volScalarField field 'interfaceEnergyFluxes'.

    All wall patches are included by default; to restrict the calculation to
    certain patches, use the optional 'patches' entry.

Usage
    Example of function object specification:
    \verbatim
    interfaceEnergyFluxes1
    {
        type        interfaceEnergyFluxes;
        libs        ("libfieldFunctionObjects.so");
        ...
        patches     (".*Wall");
    }
    \endverbatim

    Where the entries comprise:
    \table
        Property | Description                | Required   | Default value
        type     | type name: interfaceEnergyFluxes    | yes        |
        patches  | list of patches to process | no         | all wall patches
    \endtable

See also
    Foam::functionObject
    Foam::functionObjects::fvMeshFunctionObject
    Foam::functionObjects::writeFile
    Foam::functionObjects::timeControl

SourceFiles
    interfaceEnergyFluxes.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_interfaceEnergyFluxes_H
#define functionObjects_interfaceEnergyFluxes_H

#include "fvMeshFunctionObject.H"
#include "writeFile.H"
#include "volFieldsFwd.H"
#include "volFields.H"
#include "HashSet.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                       Class interfaceEnergyFluxes Declaration
\*---------------------------------------------------------------------------*/

class interfaceEnergyFluxes
:
    public fvMeshFunctionObject,
    public writeFile
{

private:

    volScalarField interfaceEnergyFluxLiquid_;
    volScalarField interfaceEnergyFluxGas_;

    scalar computeTotalEnergy
    (
        const volScalarField& energyFlux,
        const volVectorField& normal
    );

protected:
    //- File header information
    virtual void writeFileHeader(Ostream& os) const;

private:
    // Private member functions

        //- Disallow default bitwise copy construct
        interfaceEnergyFluxes(const interfaceEnergyFluxes&);

        //- Disallow default bitwise assignment
        void operator=(const interfaceEnergyFluxes&);


public:

    //- Runtime type information
    TypeName("interfaceEnergyFluxes");


    // Constructors

        //- Construct from Time and dictionary
        interfaceEnergyFluxes
        (
            const word& name,
            const Time& runTime,
            const dictionary&
        );


    //- Destructor
    virtual ~interfaceEnergyFluxes();


    // Member Functions

        //- Read the interfaceEnergyFluxes data
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
