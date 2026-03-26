/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be Ubeameful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 UbeamA


\*---------------------------------------------------------------------------*/

#include "phaseCoupling.H"
#include "mathematicalConstants.H"
#include "clock.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::phaseCoupling::checkConnection(const label startCell, const label currentID)
{
    // Set volume ID
    vofID_[startCell] = currentID;

    // Create dynamic list to store all neighbours
    DynamicList<label> cells(0);
    // Get neighbours of startCell
    cells.append(mesh_.cellCells()[startCell]);

    // Start outer loop
    while (cells.size() > 0)
    {
        // Create new dynamic list to store all new neighbours
        DynamicList<label> neighbours(0);

        // Loop over all current neighbours
        forAll(cells, i)
        {
            label cellI = cells[i];

            // Check if neighbour cell is liquid and has no ID yet
            if (alpha_[cellI] > alphaLimit_ && mag(vofID_[cellI] - currentID) > 0.1)
            {
                // Set volume ID
                vofID_[cellI] = currentID;

                // add neighbours to neighboursNeighbours list
                neighbours.append(mesh_.cellCells()[cellI]);
            }
        }
        cells.clear();
        cells.transfer(neighbours);
    }
}


void Foam::phaseCoupling::levelSetFunction()
{
    dimensionedScalar meter = dimensionedScalar("meter", dimLength, 1.0);

    // Characteristic cell length
    volScalarField dx
    (
        IOobject
        (
            "dx",
            mesh_.time().timeName(),
            mesh_
        ),
        mesh_,
        dimensionedScalar("dx", dimless, scalar(0.0)),
        zeroGradientFvPatchScalarField::typeName
    );
    dx.primitiveFieldRef() = Foam::cbrt(mesh_.V());
    dx.correctBoundaryConditions();

    // Smoothing coefficient
    dimensionedScalar nun("nun", dimArea, 0.001);

    // Coefficient for initial field
    volScalarField Gamma(0.75*dx);

    // Set initial field
    psi_ = Gamma*(2.0*alpha_ - 1.0);
    volScalarField psi0 = psi_;

    // Smoothed signed distance function
    volScalarField signPsi0 = psi0 / Foam::sqrt(sqr(psi0) + sqr(dx));

    // Artifical time step size
    volScalarField dt = 0.1*dx;

    // Number of iterations
    label nIter = 50;

    for (int iter=0; iter<nIter; iter++)
    {
        // Smoothing correction
        volVectorField n = fvc::grad(psi_) / ( mag(fvc::grad(psi_)) + SMALL/meter);
        volScalarField corr = nun*( n && fvc::grad( n &&  fvc::grad(psi_)));

        psi_ = psi_ + signPsi0*(scalar(1.0) - mag(fvc::grad(psi_)*meter) - corr)*dt;
        psi_.correctBoundaryConditions();
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


Foam::phaseCoupling::phaseCoupling
(
    const volVectorField& U,
    volScalarField& alpha,
    dropletCloud& dropletCloud
)
:
    mesh_(U.mesh()),
    dict_
    (
        IOobject
        (
            "cloudProperties",
            mesh_.time().constant(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    U_(U),
    alpha_(alpha),
    cloud_(dropletCloud),
    active_(readBool(dict_.subDict("phaseCoupling").lookup("active"))),
    alphaLimit_(readScalar(dict_.subDict("phaseCoupling").lookup("alphaLimit"))),
    dMax_(readScalar(dict_.subDict("phaseCoupling").lookup("dMax"))),
    sphMax_(readScalar(dict_.subDict("phaseCoupling").lookup("sphericity"))),
    startTime_(readScalar(dict_.subDict("phaseCoupling").lookup("startTime"))),
    nInterval_(readLabel(dict_.subDict("phaseCoupling").lookup("nInterval"))),
    vofID_
    (
        IOobject
        (
            "vofID",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("vofID", dimless, scalar(0.0)),
        zeroGradientFvPatchScalarField::typeName
    ),
    psi_
    (
        IOobject
        (
            "psi",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("psi", dimless, scalar(0.0)),
        zeroGradientFvPatchScalarField::typeName
    ),
    source_
    (
        IOobject
        (
            "source",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedVector("source", dimensionSet(0,1,-2,0,0,0,0), vector(0,0,0)),
        zeroGradientFvPatchScalarField::typeName
    ),
    damping_
    (
        IOobject
        (
            "damping",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("damping", dimensionSet(0,0,-1,0,0,0,0), scalar(0.0)),
        zeroGradientFvPatchScalarField::typeName
    )
{}

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::phaseCoupling::update()
{
    // Reset volumeID and damping/source fields
    vofID_ *= 0.0;
    damping_ *= 0.0;
    source_ *= 0.0;

    if
    (
        !active_
     || mesh_.time().timeIndex() % nInterval_ != 0
     || mesh_.time().time().value() < startTime_
    )
    {
        return;
    }


    // Global time step
    dimensionedScalar deltaT = mesh_.time().deltaT();

    // Local-to-global cell reference
    globalIndex globalNumbering (mesh_.nCells());

    // Create the volumeID field
    forAll(mesh_.cells(), cellI)
    {
        // Convert local processor cellIDs to global cellIDs
        label globalCellID = globalNumbering.toGlobal(cellI);

        // Cell is considered liquid and not yet marked
        if (alpha_[cellI] > alphaLimit_ && mag(vofID_[cellI]) < 0.1)
        {
            // Check surrounding cells
            checkConnection(cellI, globalCellID+1);
        }
    }
    vofID_.correctBoundaryConditions();


    // Correct for parallel processing
    // Number of corrections
    label corr = Foam::log(scalar(Pstream::nProcs())) / Foam::log(2.0) + 1;
    for(int i=0; i<=corr; i++)
    {
        forAll(mesh_.boundary(), patchI)
        {
            // Current patch
            const fvPatch& pPatch = mesh_.boundary()[patchI];

            // check if patch is coupled, e.g. processor patch
            if (pPatch.coupled())
            {
                // loop over all faces of processor patch
                forAll(pPatch, faceI)
                {
                    // get cell values
                    const scalar ID_own = vofID_.boundaryField()[patchI].patchInternalField()()[faceI];
                    const scalar ID_nei = vofID_.boundaryField()[patchI].patchNeighbourField()()[faceI];

                    // check, if volumeIDs on both sides are fluid
                    if ( (ID_nei > 0.1) && (ID_own > 0.1) )
                    {
                        // check, if volumeIDs on both sides are different
                        if ( mag(ID_nei - ID_own) > 0.1 )
                        {
                            scalar minID = Foam::min(ID_own, ID_nei);
                            scalar maxID = Foam::max(ID_own, ID_nei);

                            scalar diff(maxID - minID);
                            volScalarField filtered = diff
                                                  * pos( vofID_ - maxID + 0.1)
                                                  * pos(-vofID_ + maxID + 0.1);
                            vofID_ -= filtered;
                        }
                    }
                }
            }
        }
        vofID_.correctBoundaryConditions();
    }


    // Count number of continua and renumber vofID
    {
        volScalarField tmp(vofID_);
        vofID_ = 0.0;
        scalar maxID = gMax(tmp);
        int iter = 1;
        while (maxID > 0.1)
        {
            volScalarField curVolID = pos(tmp - maxID + 0.1);
            vofID_ += iter*curVolID;
            tmp -= maxID*curVolID;
            maxID = gMax(tmp);
            iter++;
        }
    }


    // Calculate droplet volume, velocity, and position
    // create lists to store data
    label maxID = floor(gMax(vofID_) + 0.5);
    labelList noCells(maxID+1, label(0));
    scalarList cellVolume(maxID+1, scalar(0));
    scalarList volume(maxID+1, scalar(0));
    vectorList position(maxID+1, vector(0,0,0));
    vectorList velocity(maxID+1, vector(0,0,0));

    // Loop over all cells and store corresponding data
    forAll(mesh_.cells(), cellI)
    {
        if (vofID_[cellI] > 0.1)
        {
            // Get volume ID
            label volID = floor(vofID_[cellI] + 0.5);

            // Store data in corresponding list
            noCells[volID] += 1;
            cellVolume[volID] += mesh_.V()[cellI];
            volume[volID] += alpha_[cellI]*mesh_.V()[cellI];
            velocity[volID] += alpha_[cellI]*mesh_.V()[cellI]*U_[cellI];
            position[volID] += alpha_[cellI]*mesh_.V()[cellI]*mesh_.C()[cellI];
        }
    }

    // Account for parallel processing
    reduce( noCells, sumOp<labelList>() );
    reduce( cellVolume, sumOp<scalarList>() );
    reduce( volume, sumOp<scalarList>() );
    reduce( position, sumOp<vectorList>() );
    reduce( velocity, sumOp<vectorList>() );

    Info << volume << endl;
    // Step F: Inject droplets
    for(int i=1; i<=maxID; i++)
    {
        if (volume[i] > VSMALL)
        {
            // Calculate droplet diameter
            const scalar d = cbrt(6.0*volume[i]/constant::mathematical::pi);

            // Droplet must be smaller than dMax_
            if ( d < dMax_ )
            {
                // Average cell size for droplet
                scalar dx = Foam::cbrt(cellVolume[i] / noCells[i]);

                // Droplet diameter has to larger than 0.5 average cell size
                // to suppress very small droplets
                if (d > 0.5*dx)
                {
                    // Weighting of position and velocity
                    vector pos = position[i] / volume[i];
                    vector U = velocity[i] / volume[i];

                    Pout << d << endl;
                    // Calculate maximal distance to center of mass for liquid structure
                    // filtered for cells only belonging to given structure
                    volScalarField mask = Foam::pos( vofID_ - i + 0.1)
                                         *Foam::pos(-vofID_ + i + 0.1);
                    vectorField cellCenters = mesh_.C();
                    scalarField distance = Foam::mag(cellCenters - pos)*mask.ref();

                    // Sphericity based on ideal sphere
                    scalar sph = 2.0*gMax(distance)/d;

                    // Sphericity limit
                    if (sph < sphMax_)
                    {
                        // Inject droplet
                        cloud_.inject(pos, d, U);

                        // Set alpha field and adjust velocity damping field
                        damping_ += GREAT*mask/deltaT;
                        alpha_ -= mask*alpha_;
                        alpha_.oldTime() -= mask*alpha_.oldTime();
                    }
                }
            }
        }
    }
}

// ************************************************************************* //
