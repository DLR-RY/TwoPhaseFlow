/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of the IsoAdvector source code library, which is an
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

Application
    isoAdvector

Description
    Advects a volume of fluid across an FVM mesh by fluxing fluid through its
    faces. Fluid transport across faces during a time step is estimated from
    the cell cutting of isosurfaces of the VOF field.

Author
    Johan Roenby, DHI, all rights reserved.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
//#include "isoAdvection.H"
#include "advectionSchemes.H"
#include "immiscibleIncompressibleTwoPhaseMixture.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createControl.H"
    #include "createTimeControls.H"

    #include "createFields.H"
    #include "readTimeControls.H"
    #include "CourantNo.H"
    #include "alphaCourantNo.H"
    #include "setDeltaT.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;
    scalar executionTime = runTime.elapsedCpuTime();

    fileName outputFile(runTime.path()/"error.dat");
    OFstream os(outputFile);

    IOdictionary isoSurfDict
    (
        IOobject
        (
            "setAlphaFieldDict",
            mesh.time().system(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    Info << "reverseTime " << reverseTime << endl;
    Info << "period " << period << endl;

    while (runTime.run())
    {
        #include "readTimeControls.H"
        #include "CourantNo.H"
        #include "alphaCourantNo.H"
        #include "setDeltaT.H"

        //Setting velocity field and face fluxes for next time step
        scalar t = runTime.time().value();
        scalar dt = runTime.deltaT().value();
        if ( reverseTime > 0.0 && t >= reverseTime )
        {
            Info<< "Reversing flow" << endl;
            phi = -phi;
            phi0 = -phi0;
            U = -U;
            U0 = -U0;
            reverseTime = -1.0;
        }
        if ( period > 0.0 )
        {
            phi = phi0*Foam::cos(2.0*M_PI*(t + 0.5*dt)/period);
            U = U0*Foam::cos(2.0*M_PI*(t + 0.5*dt)/period);
        }
        if(spirallingFlow > 0)
        {
            U = U0*Foam::cos(constant::mathematical::pi*(t+ 0.5*dt)/spirallingFlow);
            phi = phi0*Foam::cos(constant::mathematical::pi*(t+ 0.5*dt)/spirallingFlow);
        }

        runTime++;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        //Advance alpha1 from time t to t+dt
        advector->advect();

        //Write total VOF and discrepancy from original VOF to log
        label lMin = -1, lMax = -1;
        scalar aMax = -GREAT, aMin = GREAT;
        forAll(alpha1.internalField(),ci)
        {
            if ( alpha1[ci] > aMax)
            {
                aMax = alpha1[ci];
                lMax = ci;
            }
            else if ( alpha1[ci] < aMin )
            {
                aMin = alpha1[ci];
                lMin = ci;
            }
        }
        reduce(aMin, minOp<scalar>());
        reduce(aMax, maxOp<scalar>());

        const scalar V = fvc::domainIntegrate(alpha1).value();//gSum(mesh.V()*alpha1.internalField()());

        Info << "t = " << runTime.time().value() << ",\t sum(alpha*V) = " << V
             << ",\t dev = " << 100*(1.0-V/V0) << "%"
             << ",\t 1-max(alpha1) = " << 1-aMax << " at cell " << lMax
             << ",\t min(alpha1) = " << aMin << " at cell " << lMin << endl;

        runTime.write();

        //Clip and snap alpha1 to ensure strict boundedness to machine precision
      /*  if ( clipAlphaTol > 0.0 )
        {
            alpha1 = alpha1*pos(alpha1-clipAlphaTol)*neg(alpha1-(1.0-clipAlphaTol)) + pos(alpha1-(1.0-clipAlphaTol));
        }
        if ( snapAlpha )
        {
            alpha1 = min(1.0,max(0.0,alpha1));
        }*/

        scalar newExecutionTime = runTime.elapsedCpuTime();
        Info<< "ExecutionTime = " << newExecutionTime << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << "  timeStepTime = " << newExecutionTime - executionTime << " s"
            << nl << endl;
        executionTime = runTime.elapsedCpuTime();
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
