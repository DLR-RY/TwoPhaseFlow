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
    Foam::implicitFunctions::ellipsoidImplicitFunction

Description
    creates an ellipsoidImplicitFunction

    Original code supplied by Henning Scheufler, DLR (2019)

SourceFiles
    ellipsoidImplicitFunction.C

\*---------------------------------------------------------------------------*/

#ifndef implicitFunction_ellipsoidImplicitFunction_H
#define implicitFunction_ellipsoidImplicitFunction_H

#include "implicitFunction.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace implicitFunctions
{

/*---------------------------------------------------------------------------*\
                  Class ellipsoidImplicitFunction Declaration
\*---------------------------------------------------------------------------*/

class ellipsoidImplicitFunction
:
    public implicitFunction
{
private:

    // Private Member Data

        //- Axis
        vector semiAxis_;

        //- Origin point
        vector origin_;


    // Private Member Functions

        //- No copy construct
        ellipsoidImplicitFunction(const ellipsoidImplicitFunction&) = delete;

        //- No bitwise assignment
        void operator=(const ellipsoidImplicitFunction&) = delete;


public:

    //- Runtime type information
    TypeName("ellipsoidImplicitFunction");


    // Constructors

        //- Construct from components
        ellipsoidImplicitFunction(const vector& semiAxis);

        //- Construct from dictionary
        ellipsoidImplicitFunction(const dictionary& dict);


    //- Destructor
    virtual ~ellipsoidImplicitFunction() = default;

    virtual scalar value(const vector& p) const
    {
          return
           -sqrt
            (
                sqr((p.x() - origin_.x())/semiAxis_.x())
              + sqr((p.y() - origin_.y())/semiAxis_.y())
              + sqr((p.z() - origin_.z())/semiAxis_.z())
            ) + 1;
    }

    virtual vector grad(const vector& p) const
    {
        // normal_ has the length of one
        return vector
        (
            2*(p.x() - origin_.x())/sqr(semiAxis_.x()),
            2*(p.y() - origin_.y())/sqr(semiAxis_.y()),
            2*(p.z() - origin_.z())/sqr(semiAxis_.z())
        );
    }

    virtual scalar distanceToSurfaces(const vector& p) const
    {
        NotImplemented;
        return GREAT;
    }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace implicitFunctions
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
