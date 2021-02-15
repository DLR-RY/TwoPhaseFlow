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


#include "isoRDF.H"
#include "addToRunTimeSelectionTable.H"

#include "fvc.H"
#include "leastSquareInterpolate.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace reconstruction
{
    defineTypeNameAndDebug(isoRDF, 0);
    addToRunTimeSelectionTable(reconstructionSchemes,isoRDF, components);
}
}

void Foam::reconstruction::isoRDF::calcResidual
(
        const volVectorField& normal, // new
        Map<vector>& oldNormal, // old
        Map<scalar>& normalResidual
)
{

    // loop cells and calculated the gradient in that cell
    normalResidual.clear();

    forAll(normal, celli)
    {
        if (mag(normal[celli]) != 0) //
        {
            const vector cellNormal = normal[celli]/mag(normal[celli]);
            vector oldCellNormal = oldNormal[celli];
            if(oldCellNormal == vector::zero)
            {
                normalResidual.insert(celli,1);
            }
            else
            {
                oldCellNormal /= mag(oldCellNormal);
                scalar normalRes = (1 - (cellNormal & oldCellNormal));

                normalResidual.insert(celli,normalRes);
            }
        }
    }
}


void Foam::reconstruction::isoRDF::setInitPoints()
{
    interfaceLabels_.clear();
    forAll(alpha1_,celli)
    {
        if(sIterIso_.isASurfaceCell(alpha1_[celli]))
        {
            interfaceCell_[celli] = true; // is set to false earlier
            interfaceLabels_.append(celli);
        }
    }

    // resizes nextToInterfaces
    RDF_.markCellsNearSurf(interfaceCell_,1);
    const boolList& nextToInterface_ = RDF_.nextToInterface();
    exchangeFields_.updateStencil(nextToInterface_);
}


void Foam::reconstruction::isoRDF::interpolatePoints(const volScalarField& phi)
{
    leastSquareInterpolate<scalar> lsInterpol("polyDegree1",mesh_.geometricD());

    exchangeFields_.setUpCommforZone(interfaceCell_,false);

    Map<Field<vector> > mapCC(exchangePoints_.getPointDatafromOtherProc(interfaceCell_,mesh_.C()));
    Map<Field<scalar> > mapPhi(exchangePoints_.getPointDatafromOtherProc(interfaceCell_,phi));

    DynamicField<vector > cellCentre(interfaceLabels_.size());
    DynamicField<scalar > phiValues(interfaceLabels_.size());

    labelHashSet pVisisted;
    const labelListList& pCells = mesh_.cellPoints();
    const labelListList& cPoints = mesh_.pointCells();

    forAll(interfaceLabels_, i)
    {
        //if(interfaceCell_[celli])
        //{
        const label celli = interfaceLabels_[i];


        forAll(pCells[celli],j)
        {

            const label pI = pCells[celli][j];
            // point already calculated skip point
            if(pVisisted.found(pI))
            {
                continue;
            }
            pVisisted.insert(pI);
            cellCentre.clear();
            phiValues.clear();
            if(exchangePoints_.isBoundaryPoint()[pI])
            {
                cellCentre.append(mapCC[pI]);
                phiValues.append(mapPhi[pI]);
                cellCentre -= mesh_.points()[pI];
                ap_[pI] = lsInterpol.interpolate(cellCentre,phiValues);
            }
            else
            {
                forAll(cPoints[pI],k)
                {
                    const label neiCelli = cPoints[pI][k];
                    cellCentre.append(mesh_.C()[neiCelli]);
                    phiValues.append(phi[neiCelli]);
                }
                cellCentre -= mesh_.points()[pI];
                ap_[pI] = lsInterpol.interpolate(cellCentre,phiValues);
            }
        }
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::reconstruction::isoRDF::isoRDF
(
    volScalarField& alpha1,
    const surfaceScalarField& phi,
    const volVectorField& U,
    const dictionary& dict
)
:
    reconstructionSchemes
    (
        typeName,
        alpha1,
        phi,
        U,
        dict
    ),
    mesh_(alpha1.mesh()),
    ap_(mesh_.nPoints()),
    vof2IsoTol_(readScalar(modelDict().lookup("vof2IsoTol" ))),
    surfCellTol_(readScalar(modelDict().lookup("surfCellTol" ))),
    tol_(modelDict().lookupOrDefault("tol" ,1e-6)),
    iteration_(modelDict().lookupOrDefault("iterations" ,5)),
    RDF_(mesh_),
    exchangeFields_(mesh_),
    exchangePoints_(mesh_),
    sIterIso_(mesh_,ap_,surfCellTol_)


{
    Info << "iteration_ " << iteration_ << endl;
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::reconstruction::isoRDF::~isoRDF()
{}

// * * * * * * * * * * * * * * Protected Access Member Functions  * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //


// ************************************************************************* //
void Foam::reconstruction::isoRDF::reconstruct(bool forceUpdate)
{
    const bool uptodate = alreadyReconstructed(forceUpdate);

    if (uptodate && !forceUpdate)
    {
        return;
    }

    // sets interfaceLabels_ and interfaceCells_
    // update communication for RDF
    setInitPoints();
    // init ap_ with alpha field
    interpolatePoints(alpha1_);

    Map <vector> oldNormal;

    centre_ = dimensionedVector("centre",dimLength,vector::zero);
    normal_ = dimensionedVector("normal",dimArea,vector::zero);

    for(int iter=0;iter<iteration_;iter++)
    {
        forAll(interfaceLabels_, i)
        {
            const label celli = interfaceLabels_[i];
            vector n(0,0,0);
            if(mag(normal_[celli]) != 0)
            {
                n = normal_[celli]/mag(normal_[celli]);
            }

            oldNormal.set(celli,n);

            sIterIso_.vofCutCell
            (
                    celli,
                    alpha1_[celli],
                    vof2IsoTol_,
                    100
            );

            if(sIterIso_.cellStatus() == 0)
            {

                normal_[celli] = sIterIso_.surfaceArea();
                centre_[celli] = sIterIso_.surfaceCentre();
                if(mag(normal_[celli]) == 0)
                {
                    normal_[celli] = vector::zero;
                    centre_[celli] = vector::zero;
                    interfaceCell_[celli] = false;
                }


            }
            else
            {
                normal_[celli] = vector::zero;
                centre_[celli] = vector::zero;
                interfaceCell_[celli]=false;
            }
        }

        normal_.correctBoundaryConditions();
        centre_.correctBoundaryConditions();

        // nextToInterface was set in setInitNormals
        RDF_.constructRDF
        (
            RDF_.nextToInterface(),
            centre_,
            normal_,
            exchangeFields_,
            false
        );
        interpolatePoints(RDF_);
        Map<scalar> residual;
        calcResidual(normal_,oldNormal,residual);

        List<scalar> res(residual.size());
        label count = 0;
        forAllIter(Map<scalar>, residual, iter)
        {
            res[count++] = iter();
        }
        Info << "current residual absolute = " << gAverage(res) << endl;


        if(iter == 1)
        {
            Info << "intial residual absolute = " << gAverage(res) << endl;

        }

        if(( gAverage(res) < tol_ && (iter >= 1 )) || iter + 1  == iteration_)
        {
            Info << "iterations = " << iter << endl;
            Info << "final residual absolute = " << gAverage(res) << endl;
            break;
        }

    }

}
