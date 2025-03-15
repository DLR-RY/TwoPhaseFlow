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

\*---------------------------------------------------------------------------*/

#include "improvedGravity.H"
#include "fvc.H"
#include "gravityMeshObject.H"
#include "reconstructionSchemes.H"
#include "cutFaceAdvect.H"
#include "processorPolyPatch.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(improvedGravity, 0);
    addToRunTimeSelectionTable(accelerationForceModel,improvedGravity, components);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::improvedGravity::improvedGravity
(
    const dictionary& dict,
    const fvMesh& mesh
)
:
    accelerationForceModel
    (
        typeName,
        dict,
        mesh
    ),
    improvedGravityDict_(dict),
    g_
    (
        "gravity",
        dimAcceleration,
        vector(0,0,0)
    ),
    hRef_
    (
        IOobject
        (
            "hRef",
            mesh.time().constant(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        dimensionedScalar(dimLength, Zero)
    )
{
    calculateAcc();
}


// * * * * * * * * * * * * * * Public Access Member Functions  * * * * * * * //

void Foam::improvedGravity::calculateAcc()
{
    // read only if mesh changed would be clever
    const fvMesh& mesh = acc_.mesh();
    const uniformDimensionedVectorField& g = meshObjects::gravity::New(mesh.time());
    g_.value() = g.value();

    dimensionedScalar ghRef
    (
        mag(g_.value()) > SMALL
      ? g_ & (cmptMag(g_.value())/mag(g_.value()))*hRef_
      : dimensionedScalar("ghRef", g_.dimensions()*dimLength, 0)
    );

    acc_ = (g_ & mesh.C()) - ghRef;
    accf_ = (g_ & mesh.Cf()) - ghRef;


    reconstructionSchemes& surf = mesh.lookupObjectRef<reconstructionSchemes>("reconstructionScheme");
    surf.reconstruct(false);

    cutFaceAdvect advectFace(mesh, surf.alpha1());

    const volVectorField& faceCentre = surf.centre();
    const volVectorField& faceNormal = surf.normal();

    const DynamicField<label>& interfaceLabels = surf.interfaceLabels();
    boolList isSurfCell(mesh.nCells(), false);
    forAll(interfaceLabels, li)
    {
        isSurfCell[interfaceLabels[li]] = true;
    }

    const labelUList& own = mesh.owner();
    const labelUList& nei = mesh.neighbour();

    accf_ = (g & mesh.Cf()) - ghRef;

    // Setting internal accf_ values
    scalarField& accfIn = accf_.primitiveFieldRef();
    point Ci(Zero);
    forAll(accfIn, facei)
    {
        const label ownerCell(own[facei]);
        const label neiCell(nei[facei]);
        if (isSurfCell[ownerCell])
        {
            Ci = faceCentre[ownerCell];
            if (isSurfCell[neiCell])
            {
                Ci += faceCentre[neiCell];
                Ci *= 0.5;
            }
        }
        else if (isSurfCell[neiCell])
        {
            Ci = faceCentre[neiCell];
        }
        accfIn[facei] = (g.value() & Ci) -  ghRef.value();
    }

    // Setting boundary accf_ values
    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    forAll(patches, patchi)
    {
        // Pout << "patches[patchi].name() = " << patches[patchi].name() << endl;
        scalarField& accfbi = accf_.boundaryFieldRef()[patchi];
        forAll(accfbi, facei)
        {
            const label globalFacei = facei + patches[patchi].start();
            const label ownerCell(own[globalFacei]);

            if (isSurfCell[ownerCell])
            {
                // Reconstructing height on patch from intersection with interface
                // in owner cell
                advectFace.calcSubFace
                (
                    globalFacei,
                    faceNormal[ownerCell],
                    faceCentre[ownerCell]
                );
                DynamicList<point> edgePoints(advectFace.surfacePoints());
                point Ci(Zero);
                forAll(edgePoints, pointi)
                {
                    Ci += edgePoints[pointi];
                }
                if (edgePoints.size() == 2)
                {
                    Ci /= edgePoints.size();
                    Info << "Info: face " << globalFacei << " has two edge points" << endl;
                    Info << "Info: owner cell " << ownerCell << endl;
                    Info << "Info: Ci " << Ci << endl;
                }
                else
                {
                    Ci = faceCentre[ownerCell];
                    Info << "Warning: Reconstructed face-interface intersection line on face " << globalFacei << " has points: " << edgePoints << endl;
                    Info << "Info: face " << globalFacei << " has two edge points" << endl;
                    Info << "Info: owner cell " << ownerCell << endl;
                    Info << "Info: Ci " << Ci << endl;
                }
                accfbi[facei] = (g.value() & Ci) - ghRef.value();
            }
            //  else special case where owner has alpha = 1 and adjacent boundary
            // cell has alpha = 0. Use Cf of shared face but how?
        }
    }

    // Correcting ghf on cyclic patches
    // Loop over all patches
    forAll(patches, patchi)
    {
        const fvPatch& patch = mesh.boundary()[patchi];

        // Check if the patch is cyclic
        if (isA<cyclicFvPatch>(patch))
        {
    //        Info << "Patch " << patchi << " is cyclic" << endl;
            const cyclicFvPatch& cyclicPatch = refCast<const cyclicFvPatch>(patch);

            // Access field values on the patches
            scalarField& patchField = accf_.boundaryFieldRef()[patchi];
            scalarField& neighbPatchField = accf_.boundaryFieldRef()[cyclicPatch.neighbPatchID()];

            // Set field values to the mean of the corresponding face pairs
            forAll(patchField, facei)
            {
                scalar faceValue = 0.5 * (patchField[facei] + neighbPatchField[facei]);
                patchField[facei] = faceValue;
                neighbPatchField[facei] = faceValue;
            }
        }
    }

    // Setting accf_ on processor patches
    if (Pstream::parRun())
    {
        DynamicList<label> neighProcs;
        PstreamBuffers pBufs(Pstream::commsTypes::nonBlocking);

        // Send
        for (const polyPatch& pp : patches)
        {
            const auto* ppp = isA<processorPolyPatch>(pp);

            if (ppp && pp.nPoints())
            {
                const auto& procPatch = *ppp;
                const label nbrProci = procPatch.neighbProcNo();

                neighProcs.append(nbrProci);
                UOPstream toNbr(nbrProci, pBufs);


                const scalarField& accfbi = accf_.boundaryField()[procPatch.index()];
                List<label> interfaceFaces;
                List<scalar> interfaceFaceValues;
                forAll(accfbi, facei)
                {
                    const label globalFacei = facei + procPatch.start();
                    const label ownerCell(own[globalFacei]);
                    if (isSurfCell[ownerCell])
                    {
                        interfaceFaces.append(facei);
                        interfaceFaceValues.append(accfbi[facei]);
                    }
                }

                /*
                Pout << "Processor " << Pstream::myProcNo() <<
                    " sending interfaceFaces = " << interfaceFaces <<
                    " and interfaceFaceValues = " << interfaceFaceValues <<
                    " to neighbour processor = " << nbrProci << endl;
                */

                toNbr << interfaceFaces << interfaceFaceValues;
            }
        }

        // Limited to involved neighbour procs
        pBufs.finishedNeighbourSends(neighProcs);


        // Receive and combine
        for (const polyPatch& pp : patches)
        {
            const auto* ppp = isA<processorPolyPatch>(pp);

            if (ppp && pp.nPoints())
            {
                const auto& procPatch = *ppp;

                const label nbrProci = procPatch.neighbProcNo();

                UIPstream fromNeighb(nbrProci, pBufs);
                List<label> neighbourInterfaceFaces;
                List<scalar> neighbourInterfaceFaceValues;

                fromNeighb >> neighbourInterfaceFaces >> neighbourInterfaceFaceValues;

                /*
                Pout << "Processor " << Pstream::myProcNo() <<
                    " received  neighbourInterfaceFaces = " <<  neighbourInterfaceFaces <<
                    " and neighbourInterfaceFaceValues = " << neighbourInterfaceFaceValues <<
                    " from neighbour processor: " << nbrProci << endl;
                */

                scalarField& accfbi = accf_.boundaryFieldRef()[procPatch.index()];

                forAll(neighbourInterfaceFaces, fi)
                {
                    const label facei = neighbourInterfaceFaces[fi];
                    accfbi[facei] += neighbourInterfaceFaceValues[fi];
                    const label globalFacei = facei + procPatch.start();
                    const label ownerCell(own[globalFacei]);
                    if (isSurfCell[ownerCell])
                    {
                        accfbi[facei] *= 0.5;
                    }
                }
            }
        }
    }
}


Foam::tmp<Foam::surfaceScalarField> Foam::improvedGravity::accelerationForce()
{
    const fvMesh& mesh = acc_.mesh();
    const volScalarField& rho = mesh.lookupObject<volScalarField>("rho");
    return -accf()*fvc::snGrad(rho);
}

// ************************************************************************* //
