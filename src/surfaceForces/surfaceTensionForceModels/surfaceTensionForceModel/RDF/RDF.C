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

#include "RDF.H"
#include "addToRunTimeSelectionTable.H"

#include "alphaContactAngleTwoPhaseFvPatchScalarField.H"
#include "mathematicalConstants.H"
#include "surfaceInterpolate.H"
#include "fvcDiv.H"
#include "fvcGrad.H"
#include "fvcSnGrad.H"
#include "fvc.H"

#include "plane.H"
#include "interpolationCellPoint.H"

#include "reconstructionSchemes.H"
#include "wedgePolyPatch.H"
#include "indexedOctree.H"
#include "treeDataPoint.H"
#include "processorPolyPatch.H"
#include "syncTools.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(RDF, 0);
    addToRunTimeSelectionTable(surfaceTensionForceModel,RDF, components);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::RDF::RDF
(
    const dictionary& dict,
    const volScalarField& alpha1,
    const surfaceScalarField& phi,
    const volVectorField& U
)
:
    surfaceTensionForceModel
    (
        typeName,
        dict,
        alpha1,
        phi,
        U
    ),
    curvFromTr_(dict.lookupOrDefault("curvFromTr",true)),
    deltaN_
    (
        "deltaN",
        1e-8/pow(average(alpha1.mesh().V()), 1.0/3.0)
    ),

    FaceCentreField_(alpha1.mesh().nCells()),
    FaceNormalField_(alpha1.mesh().nCells()),
    RDF_(reconstructedDistanceFunction::New(alpha1.mesh()))
{
    // curvature from trace does not work with wedges
    bool wedge = false;
    const polyBoundaryMesh& patches = alpha1.mesh().boundaryMesh();
    for (const polyPatch& pp : patches)
    {
        if(isA<wedgePolyPatch>(pp))
        {
            wedge = true;
        }
    }

    if(wedge)
    {
        curvFromTr_ = false;
    }
}

// * * * * * * * * * * * * * * Public Access Member Functions  * * * * * * * //
template <typename T>
void Foam::RDF::distributeField(const labelList& neiProcs,Field<T>& field)
{
    Field<Field <T>> sendFields(Pstream::nProcs());

    forAll(neiProcs,i)
    {
        label neiProc = neiProcs[i];
        sendFields[neiProc] = field;
    }

    Field<Field <T>> recvFields(Pstream::nProcs());

    Pstream::exchange<Field<T>, T>(sendFields, recvFields);

    forAll(recvFields,i)
    {
        field.append(recvFields[i]);
    }

}


Foam::labelList Foam::RDF::getNeibourProcs(const boolList& nextToInterface)
{
    const polyBoundaryMesh& patches = alpha1_.mesh().boundaryMesh();
    const fvMesh& mesh = alpha1_.mesh();
    const labelListList& cPoints = mesh.pointCells();
    Map<labelList> neiProcMap;
    labelList myProc(1,Pstream::myProcNo());

    for (const polyPatch& pp : patches)
    {
        if (isA<processorPolyPatch>(pp))
        {
            const labelList& pPoints = pp.meshPoints();

            for (label pI: pPoints)
            {
                bool interfacePoint = false;
                forAll(cPoints[pI],i)
                {
                    if (nextToInterface[cPoints[pI][i]])
                    {
                        interfacePoint = true;
                        break;
                    }
                }
                if (interfacePoint)
                {
                    neiProcMap.insert(pI,myProc);
                }
            }
        }
    }

    syncTools::syncPointMap(alpha1_.mesh(), neiProcMap, appendOp<label>());
    labelHashSet neiProcs;

    forAllConstIters(neiProcMap,iter)
    {
        neiProcs.insert(iter());
    }

    return neiProcs.toc();
}


Foam::label Foam::RDF::closestDistToSurface(const point& p)
{
    const label& minLabel = findMin(mag(FaceCentreField_-p)());

    return minLabel;

}


Foam::scalar Foam::RDF::distanceToSurfacePlane(const point& p)
{
    const label& minLabel = findMin(mag(FaceCentreField_-p)());

    vector c1 = FaceCentreField_[minLabel];


    vector n1 = FaceNormalField_[minLabel];

   plane nearSurface(c1,n1);
   scalar dist=mag(nearSurface.distance(p));


      if (((c1-p) & n1) > 0)
      {
         dist*=-1;
      }


  return dist;

}

void Foam::RDF::correctContactAngle
(
    surfaceVectorField::Boundary& nHatb,
    surfaceVectorField::Boundary& gradAlphaf
)
{
    scalar convertToRad = Foam::constant::mathematical::pi/180.0;

    const fvMesh& mesh = alpha1_.mesh();
    const volScalarField::Boundary& abf = alpha1_.boundaryField();

    const fvBoundaryMesh& boundary = mesh.boundary();

    forAll(boundary, patchi)
    {
        if (isA<alphaContactAngleTwoPhaseFvPatchScalarField>(abf[patchi]))
        {
            alphaContactAngleTwoPhaseFvPatchScalarField& acap =
                const_cast<alphaContactAngleTwoPhaseFvPatchScalarField&>
                (
                    refCast<const alphaContactAngleTwoPhaseFvPatchScalarField>
                    (
                        abf[patchi]
                    )
                );

            fvsPatchVectorField& nHatp = nHatb[patchi];
            const scalarField theta
            (
                convertToRad*acap.theta(U_.boundaryField()[patchi], nHatp)
            );

            const vectorField nf
            (
                boundary[patchi].nf()
            );

            // Reset nHatp to correspond to the contact angle

            const scalarField a12(nHatp & nf);
            const scalarField b1(cos(theta));

            scalarField b2(nHatp.size());
            forAll(b2, facei)
            {
                b2[facei] = cos(acos(a12[facei]) - theta[facei]);
            }

            const scalarField det(1.0 - a12*a12);

            scalarField a((b1 - a12*b2)/det);
            scalarField b((b2 - a12*b1)/det);

            nHatp = a*nf + b*nHatp;
            nHatp /= (mag(nHatp) + deltaN_.value());

            acap.gradient() = (nf & nHatp)*mag(gradAlphaf[patchi]);
            acap.evaluate();
        }
    }
}


void Foam::RDF::correct()
{
    deltaFunctionModel_->correct();

    const fvMesh& mesh = alpha1_.mesh();
    mesh.time().cpuTimeIncrement();
    const surfaceVectorField& Sf = mesh.Sf();

    reconstructionSchemes& surf =
        mesh.lookupObjectRef<reconstructionSchemes>("reconstructionScheme");

    surf.reconstruct(false);

    const volVectorField& faceCentre = surf.centre();
    const volVectorField& faceNormal = surf.normal();
    boolList interfaceCells = surf.interfaceCell();


    FaceCentreField_.clear();
    FaceNormalField_.clear();

    FaceCentreField_.setCapacity(mesh.nCells());
    FaceNormalField_.setCapacity(mesh.nCells());
    DynamicField<scalar> KatInterFace(mesh.nCells());

    forAll(faceCentre,cellI)
    {
        if (mag(faceNormal[cellI]) != 0)
        {
            FaceCentreField_.append(faceCentre[cellI]);
            FaceNormalField_.append(faceNormal[cellI]);
            interfaceCells[cellI] = true;
        }
        else
        {
            interfaceCells[cellI] = false;
        }
    }

    RDF_.markCellsNearSurf(interfaceCells,3);

    labelList neiProcs = getNeibourProcs(RDF_.nextToInterface());

    distributeField<vector>(neiProcs,FaceCentreField_);
    distributeField<vector>(neiProcs,FaceNormalField_);

    RDF_.constructRDFOctree(RDF_.nextToInterface(),FaceCentreField_,FaceNormalField_);

    RDF_.correctBoundaryConditions(); // update proc patches

    dimensionedScalar deltaN
    (
        "deltaN",
        1e-14/pow(average(mesh.V()), 1.0/3.0)
    );

    volVectorField gradRDF(fvc::grad(RDF_));

    gradRDF /= (mag(gradRDF)+deltaN*dimensionedScalar("0", dimLength, 1));

    gradRDF.correctBoundaryConditions();

    surfaceVectorField interfaceVec("interfaceVec",fvc::interpolate((gradRDF)));

    surfaceVectorField normalVec("normalVec",interfaceVec/(mag(interfaceVec)+deltaN*dimensionedScalar("0", dimLength, 1))); // macht vielleicht den fehler

    // correct contact angle
    correctContactAngle(normalVec.boundaryFieldRef(), interfaceVec.boundaryFieldRef());


    if (curvFromTr_)
    {
        const fvBoundaryMesh& boundary = mesh.boundary();

        forAll(boundary, patchi)
        {
                fvPatchVectorField& nHatp = gradRDF.boundaryFieldRef()[patchi];
                nHatp = normalVec.boundaryFieldRef()[patchi];
        }
    }

    // Face unit interface normal flux
    nHatf_ = normalVec & Sf;

    // Simple expression for curvature
    if (curvFromTr_)
    {
        K_ = -tr(fvc::grad(gradRDF));
    }
    else
    {
        K_ = -fvc::div(nHatf_);
    }

   volScalarField curvature("curvature",K_);

   interpolationCellPoint<scalar> interpolCurv(curvature);

    forAll(K_,cellI)
    {
           if (mag(faceNormal[cellI]) != 0)
           {
               K_[cellI]=  interpolCurv.interpolate(faceCentre[cellI],cellI,-1);
               KatInterFace.append(K_[cellI]);
           }
           else
           {
               K_[cellI]= 0;
           }
    }

    distributeField<scalar>(neiProcs,KatInterFace);

    Random rndGen(17301893);

    // Slightly extended bb. Slightly off-centred just so on symmetric
    // geometry there are less face/edge aligned items.
    treeBoundBox bb
    (
        treeBoundBox(FaceCentreField_).extend(rndGen, 1e-4)
    );

    bb.min() -= point(1e-8, 1e-8, 1e-8);
    bb.max() += point(1e-8, 1e-8, 1e-8);

    indexedOctree<treeDataPoint> surfaceTree
    (
        treeDataPoint
        (
            FaceCentreField_
        ),
        bb,     // bb
        8,      // maxLevel
        10,     // leafsize
        3.0     // duplicity
    );

    forAll(K_,cellI)
    {
        if (mag(faceNormal[cellI]) == 0 && RDF_.nextToInterface()[cellI])
        {
            pointIndexHit pHit =  surfaceTree.findNearest(mesh.C()[cellI], GREAT);
            const label idx = pHit.index();
            if (idx != -1)
            {
                K_[cellI]=  KatInterFace[idx];
            }
        }
        if (!RDF_.nextToInterface()[cellI])
        {
            K_[cellI]=  0;
        }
    }

    K_.correctBoundaryConditions();

    Kf_ = fvc::interpolate(K_);
}



// ************************************************************************* //
