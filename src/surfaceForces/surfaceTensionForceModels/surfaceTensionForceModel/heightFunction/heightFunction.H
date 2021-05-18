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
    Foam::heightFunction

Description
    Height Function method for the curvature compuation only works on grids
    with cuboid cells

SourceFiles
    heightFunction.C

\*---------------------------------------------------------------------------*/

#ifndef heightFunction_H
#define heightFunction_H

#include "typeInfo.H"
#include "volFields.H"
#include "dimensionedScalar.H"
#include "autoPtr.H"

#include "surfaceTensionForceModel.H"
#include "HFStencil.H"
#include "markInterfaceRegion.H"
#include "twoDimFDStencil.H"
#include "HashTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class heightFunction Declaration
\*---------------------------------------------------------------------------*/

class heightFunction
    : public surfaceTensionForceModel
{

private:

    //-
    //- Stabilisation for normalisation of the interface normal
    typedef HashTable<twoDimFDStencil,Vector2D<label>,Hash<Vector2D<label>>> twoDimStencilMap;

    const dimensionedScalar deltaN_;

    const fvMesh& mesh_;

    HFStencil stencilHF_;

    markInterfaceRegion IFRegion_;

    bool twoDim_;

    bool fullColumn(const scalar avgHeight,const scalar tol);

    void nextCell
    (
        const label celli,
        const HFStencil::orientation orientation,
        twoDimFDStencil& HFCol
    );

    void exchangeStencils
    (
        const List<List<twoDimFDStencil>>& sendStencil,
        List<twoDimFDStencil>& recvStencil
    );

    // Private Member Functions

        virtual void correctContactAngle
        (
            surfaceVectorField::Boundary& nHat,
            surfaceVectorField::Boundary& heightFunctionf
        );

        void computeColumns
        (
            const label dir,
            const Map<scalar>& mapAlpha,
            const label celli,
            const HFStencil::orientation orientation,
            twoDimFDStencil& HFCol
        );

        void getStencilValues
        (
            const Map<scalar>& MapAlpha,
            const labelList& stencil,
            DynamicField<scalar>& alphaValues
        );

        scalar calcCurvature
        (
            const scalarField& fit
        );

        //- Disallow copy construct
        heightFunction(const heightFunction&);

        //- Disallow default bitwise assignment
        void operator=(const heightFunction&);

public:

    //- Runtime type information
    TypeName("heightFunction");


    // Constructors

        //- Construct from components
        heightFunction
        (
            const dictionary& dict,
            const volScalarField& alpha1,
            const surfaceScalarField& phi,
            const volVectorField& U
        );

    //- Destructor
    virtual ~heightFunction()
    {}

    virtual void correct();


};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
