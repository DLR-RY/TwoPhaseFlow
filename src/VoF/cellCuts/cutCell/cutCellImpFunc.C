/*---------------------------------------------------------------------------*\
|             isoAdvector | Copyright (C) 2016 Johan Roenby, DHI              |
-------------------------------------------------------------------------------

License
    This file is part of IsoAdvector, which is an unofficial extension to
    OpenFOAM.

    IsoAdvector is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    IsoAdvector is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with IsoAdvector. If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "cutCellImpFunc.H"
#include "triSurfaceTools.H"
#include "SortableList.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


Foam::point Foam::cutCellImpFunc::bisection(point p0, point p1)
{
    // return
    //vector start = p0;
    //vector end = p1;
    scalar a = func_.value(p0);
    scalar b = func_.value(p1);
    if (a * b > 0)
    {
        return vector(GREAT,GREAT,GREAT);
    }

    if(a == 0)
    {
        return p0;
    }

    if(b == 0)
    {
        return p1;
    }

    scalar c = a;
    vector p3 = p0;
    while (mag(p1-p0) >= SMALL)
    {
        // Find middle point
        p3 = 0.5*(p0+p1);
        c = func_.value(p3);

        // Check if middle point is root
        if (c == 0.0)
        {
            break;


        // Decide the side to repeat the steps
        }
        else if (c*a < 0)
        {
            b = c;
            p1 = p3;
        }
        else
        {
            a = c;
            p0 = p3;
        }
    }

    return p3;
}

void Foam::cutCellImpFunc::getBase(const vector& n, vector& e0, vector& e1)
const
{
    // Guess for vector normal to n.
    vector base(1,0,0);

    scalar nComp = n & base;

    if (mag(nComp) > 0.8)
    {
        // Was bad guess. Try with different vector.

        base.x() = 0;
        base.y() = 1;

        nComp = n & base;

        if (mag(nComp) > 0.8)
        {
            base.y() = 0;
            base.z() = 1;

            nComp = n & base;
        }
    }


    // Use component normal to n as base vector.
    e0 = normalised(base - nComp * n);

    e1 = n ^ e0;

    //Pout<< "Coord system:" << endl
    //    << "    n  : " << n << ' ' << mag(n) << endl
    //    << "    e0 : " << e0 << ' ' << mag(e0) << endl
    //    << "    e1 : " << e1 << ' ' << mag(e1) << endl
    //    << endl;
}

Foam::face Foam::cutCellImpFunc::sortInterfacePoints(const DynamicList< point>& interfacePoints,vector n)
{
    vector e0, e1;
    n /= mag(n);
    getBase(n, e0, e1);

    point ctr = average(interfacePoints);
    // Get sorted angles from point on loop to centre of loop.
    SortableList<scalar> sortedAngles(interfacePoints.size());

    forAll(sortedAngles, i)
    {
        vector toCtr = normalised(interfacePoints[i] - ctr);
        tensor projPlane = tensor::I - n*n;
        toCtr = toCtr & projPlane;

        sortedAngles[i] = pseudoAngle(e0, e1, toCtr);
    }
    sortedAngles.sort();
    // Info << sortedAngles.indices() << endl;

    face f(sortedAngles.indices());

    return f;
}
Foam::point Foam::cutCellImpFunc::newtonMethod(const point& p0)
{
    vector p1 = p0;
    scalar funcValue = func_.value(p1);
    label iter = 0;
    while(mag(funcValue) >=SMALL && iter < 100)
    {
        vector gradFunc = func_.grad(p1);
        funcValue = func_.value(p1);
        p1 -= funcValue*gradFunc/mag(gradFunc);
        iter++;
    }


    return p1;
}

Foam::point Foam::cutCellImpFunc::newtonMethod(const tensor& t,const point& p0)
{
    vector p1 = p0;
    scalar funcValue = func_.value(p1);
    label iter = 0;
    while(mag(funcValue) >=SMALL && iter < 100)
    {
        vector gradFunc = func_.grad(p1);
        gradFunc = t & gradFunc ;
        if(mag(gradFunc) == 0)
        {
            break; // needs better error handling
        }
        funcValue = func_.value(p1);
        p1 -= funcValue*gradFunc/mag(gradFunc);
        iter++;
    }

    return p1;
}

Foam::triSurface Foam::cutCellImpFunc::refineFreeSurfaces(triSurface& triF)
{

    // for(int counter = 0;counter<3;counter++)
    {
        label oldSize = triF.points().size();
        triF = triSurfaceTools::redGreenRefine
        (
            triF,
            identity(triF.size())  //Hack: refine all
        );
        pointField newPoints(triF.points().size());

        for(label pI = 0; pI<triF.points().size() ;pI++)
        {
            if(pI >= oldSize)
            {
                newPoints[pI] = newtonMethod(triF.points()[pI]);
            }
            else
            {
                newPoints[pI] = triF.points()[pI];
            }
        }

        triF.movePoints(newPoints);

//        for(label pI = 0; pI<triF.points().size() ;pI++)
//        {
//            ioc.findInside(newPoints[pI]);
//        }

    }
    return triF;
}

Foam::DynamicList< Foam::point > Foam::cutCellImpFunc::refineCuttedFace
(
    const vector& n,
    DynamicList< point>& cuttedFace,
    DynamicList< point>& interfacePoints
)
{
    DynamicList< point >refinedFace(cuttedFace.size()+3);
    interfacePoints.clear();

    forAll(cuttedFace,i)
    {
        const label nextIdx = (i+1) % cuttedFace.size();
        scalar fVnextIdx = mag(func_.value(cuttedFace[nextIdx]));
        scalar fVIdx = mag(func_.value(cuttedFace[i]));

        // we found a surface edge
        refinedFace.append(cuttedFace[i]);
        if(fVIdx <= 1e-10 && fVnextIdx <= 1e-10)
        {
            point refPoint = 0.5*cuttedFace[i] + 0.5*cuttedFace[nextIdx];
            tensor t = tensor::I - n*n;
            refPoint = newtonMethod(t,refPoint);
            refinedFace.append(refPoint);
        }

    }
    forAll(refinedFace,i)
    {
        scalar fV = mag(func_.value(refinedFace[i]));
        if(fV <= 1e-10)
        {
            interfacePoints.append(refinedFace[i]);
        }
    }

    return refinedFace;
}


Foam::cutCellImpFunc::cutCellImpFunc(const fvMesh& mesh, scalarField& f, implicitFunction& func)
:
      cutCell(mesh),
      mesh_(mesh),
      cellI_(-1),
      f_(f),
      func_(func),
      cutValue_(0),
      cutFace_(cutFaceImpFunc(mesh_, f_,func)),
      cutFaceCentres_(10),
      cutFaceAreas_(10),
      isoFaceEdges_(10),
      facePoints_(10),
      faceCentre_(vector::zero),
      faceArea_(vector::zero),
      subCellCentre_(vector::zero),
      subCellVolume_(-10),
      VOF_(-10),
      cellStatus_(-1),
      freeSurf_()
{
    clearStorage();
    forAll(f_,pI)
    {
        f_[pI] = func_.value(mesh_.points()[pI]);
    }
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::label Foam::cutCellImpFunc::calcSubCell
(
    const label& cellI,
    const scalar cutValue
)
{

    // resets data members
    clearStorage();
    cellI_ = cellI;
    cutValue_ = cutValue;
    const cell& c = mesh_.cells()[cellI_];

    bool fullyBelow = true;
    bool fullyAbove = true;
    scalar avgCutArea = 0;
    List<DynamicList < point> > cutFacePoints(c.size());
    DynamicList< label > fluidFace(c.size());
    List<  vector > cutFaceNormal(c.size());

    // loop over cell faces
    forAll(c, fi)
    {
        const label facei = c[fi];

        const label faceStatus = cutFace_.calcSubFace(facei, cutValue_);

        if (faceStatus == 0) // face is cut
        {
            cutFaceCentres_.append(cutFace_.subFaceCentre());
            cutFaceAreas_.append(cutFace_.subFaceArea());
            isoFaceEdges_.append(cutFace_.surfacePoints());
            avgCutArea += mag(cutFace_.subFaceArea())/mag(mesh_.faceAreas()[facei]);
            fullyBelow = false;
            fullyAbove = false;
            cutFacePoints[fi] = cutFace_.subFacePoints();
            cutFaceNormal[fi] = cutFace_.subFaceArea();
        }
        else if (faceStatus == -1) // face fully below
        {
            cutFaceCentres_.append(cutFace_.subFaceCentre());
            cutFaceAreas_.append(cutFace_.subFaceArea());
            avgCutArea += mag(cutFace_.subFaceArea())/mag(mesh_.faceAreas()[facei]);
            fluidFace.append(facei);
            fullyAbove = false;
        }
        else
        {
            fullyBelow = false;
        }
    }
    avgCutArea /= c.size();

    if (!fullyBelow && !fullyAbove) // cell cut at least at one face
    {
        cellStatus_ = 0;


        // calc faceArea and faceCentre
        calcGeomDataCutFace
        (
            isoFaceEdges_,
            average(cutFaceCentres_),
            faceArea_,
            faceCentre_
        );

        if(mag(faceArea_) == 0)
        {

            if(avgCutArea > 0.5)
            {
                cellStatus_ = -1;
                subCellCentre_ = mesh_.C()[cellI_];
                subCellVolume_ = mesh_.V()[cellI_];
                VOF_ = 1;
            }
            else
            {
                cellStatus_ = 1;
                subCellCentre_ = vector::zero;
                subCellVolume_ = 0;
                VOF_ = 0;
            }
            return cellStatus_;
        }

        // pointField surfPoints(facePoints());
        // face freeF (identity(surfPoints.size()));
        // DynamicList<face> fL;

        // freeF.triangles (surfPoints, fL);

        // triFaceList triF(fL.size());
        // forAll(triF,i)
        // {
        //     triF[i] = triFace(fL[i][0],fL[i][1],fL[i][2]);
        // }

        // freeSurf_ = triSurface(triF,surfPoints);

        // forAll(freeSurf_,triI)
        // {
        //     cutFaceCentres_.append(freeSurf_.Sf()[triI]);
        //     cutFaceAreas_.append(freeSurf_.Cf()[triI]);
        // }

        // calc volume and sub cell centre
        calcCellData
        (
            cutFaceCentres_,
            cutFaceAreas_,
            subCellCentre_,
            subCellVolume_
        );

        // switch orientation that it points in the gas phase
        if ((faceArea_ & (faceCentre_ - subCellCentre_)) < 0)
        {
            faceArea_ *= (-1);
        }

        VOF_ = subCellVolume_ / mesh_.V()[cellI_];

        // DEBUGGING
        //scalar error = GREAT;
        
        DynamicList<point> interfacePoints(1000);
        DynamicList<point> interfaceLocalFacePoints(1000);


        for(label refineCounter = 0; refineCounter < 5; refineCounter++)
        {
            // freeSurf_ = refineFreeSurfaces(freeSurf_);

            // fileName fileFreeSurf = "freeSurf" +  std::to_string(refineCounter) +  ".vtk";
            // freeSurf_.write(fileFreeSurf);
            cutFaceCentres_.clear();
            cutFaceAreas_.clear();
            interfacePoints.clear();
            forAll(fluidFace,ffI)
            {
                const label faceI = fluidFace[ffI];
                cutFaceCentres_.append(mesh_.faceCentres()[faceI]);
                cutFaceAreas_.append(mesh_.faceAreas()[faceI]);
            }

            forAll(cutFacePoints,i)
            {
                if(cutFacePoints[i].size() != 0)
                {
                    // pointField cfP(cutFacePoints[i]);
                    vector n = cutFaceNormal[i];
                    if(mag(n) != 0)
                    {
                        n /= mag(n);
                    }
                    pointField cfP = pointField(refineCuttedFace(n,cutFacePoints[i],interfaceLocalFacePoints));
                    interfacePoints.append(interfaceLocalFacePoints);
                    face cuttedFace (identity(cfP.size()));
                    cutFacePoints[i] = DynamicList <point>(cfP);
                    label nPoints = cutFacePoints[i].size();
                    if(nPoints == 3)
                    {
                        cutFaceCentres_.append(cuttedFace.centre(cfP));
                        cutFaceAreas_.append(cuttedFace.areaNormal(cfP));
                    }
                    else
                    {
                        forAll(cutFacePoints[i],idx)
                        {
                            vector centre = cuttedFace.centre(cfP);
                            const label nextIdx = (idx+1) % cutFacePoints[i].size();
                            triPointRef triangle
                            (
                                cutFacePoints[i][idx],
                                cutFacePoints[i][nextIdx],
                                centre
                            );

                            cutFaceCentres_.append(triangle.centre());
                            cutFaceAreas_.append(triangle.areaNormal());

                        }
                    }


                    // postProcessing
                    //face cuttedFace (identity(cfP.size()));
                    if(refineCounter == 5-1)
                    {
                        DynamicList<face> fL;

                        cuttedFace.triangles (cfP, fL);

                        triFaceList triF(fL.size());
                        forAll(triF,i)
                        {
                            triF[i] = triFace(fL[i][0],fL[i][1],fL[i][2]);
                        }
                        triSurface test(triF,cfP);
                        fileName file = "cutFaceCounter_cellI" + std::to_string(cellI) + "i" + std::to_string(i) +  +  ".vtk";
                        // test.write(file);
                    }
                }
            }

        }

            // sort interface points

        point centre = average(interfacePoints);

        DynamicList<face> fLIf;
        face IFFace = sortInterfacePoints(interfacePoints,faceArea_);
        //face IFFace (identity(interfacePoints.size()));
        centre= newtonMethod(centre);
        interfacePoints.append(centre);

        //IFFace.triangles (interfacePoints, fLIf);

        triFaceList triFIF(IFFace.size());
        label centreIdx = interfacePoints.size() -1;
        forAll(triFIF,i)
        {
            const label nextIdx = (i+1) % IFFace.size();
            triFIF[i] = triFace(centreIdx,IFFace[i],IFFace[nextIdx]);
        }
        pointField IFPoints(interfacePoints);
        freeSurf_ = triSurface(triFIF,IFPoints);
        fileName file = "triiFFace_celli"  + std::to_string(cellI) +  +  ".vtk";
        for(label refineCounter = 0; refineCounter < 3; refineCounter++)
        {
            freeSurf_ = refineFreeSurfaces(freeSurf_);
        }
        // freeSurf_.write(file);

        forAll(freeSurf_,triI)
        {
            cutFaceCentres_.append(freeSurf_.Cf()[triI]);
            cutFaceAreas_.append(freeSurf_.Sf()[triI]);
        }

        // calc volume and sub cell centre
        calcCellData
        (
            cutFaceCentres_,
            cutFaceAreas_,
            subCellCentre_,
            subCellVolume_
        );
        scalar newVoF = subCellVolume_ / mesh_.V()[cellI_];
        VOF_ = newVoF;

        // DEBUGGING
        // error = (newVoF-VOF_)/newVoF;
        // Info << "error " << error << endl;
        // Info << "celli " << cellI <<  " VOF_ " << VOF_ << endl;


    }
    else if (fullyAbove) // cell fully above isosurface
    {
        cellStatus_ = 1;
        subCellCentre_ = vector::zero;
        subCellVolume_ = 0;
        VOF_ = 0;
    }
    else if (fullyBelow) // cell fully below isosurface
    {
        cellStatus_ = -1;
        subCellCentre_ = mesh_.C()[cellI_];
        subCellVolume_ = mesh_.V()[cellI_];
        VOF_ = 1;
    }


    return cellStatus_;
}

Foam::point Foam::cutCellImpFunc::subCellCentre()
{
    return subCellCentre_;
}

Foam::scalar Foam::cutCellImpFunc::subCellVolume()
{
    return subCellVolume_;
}

Foam::DynamicList<Foam::point> Foam::cutCellImpFunc::facePoints()
{
    if (facePoints_.size() == 0)
    {
        // get face points in sorted order
        calcIsoFacePointsFromEdges
        (
            faceArea_,
            faceCentre_,
            isoFaceEdges_,
            facePoints_
        );
    }

    return facePoints_;
}

Foam::point Foam::cutCellImpFunc::faceCentre()
{
    return faceCentre_;
}

Foam::vector Foam::cutCellImpFunc::faceArea()
{
    return faceArea_;
}
Foam::label Foam::cutCellImpFunc::cellStatus()
{
    return cellStatus_;
}

Foam::scalar Foam::cutCellImpFunc::VolumeOfFluid()
{
    return VOF_;
}

Foam::scalar Foam::cutCellImpFunc::cutValue() const
{
    return cutValue_;
}

const Foam::triSurface& Foam::cutCellImpFunc::triFreeSurf()
{
    return freeSurf_;
}

void Foam::cutCellImpFunc::clearStorage()
{
    cellI_ = -1;
    cutValue_ = 0;
    cutFaceCentres_.clear();
    cutFaceAreas_.clear();
    isoFaceEdges_.clear();
    facePoints_.clear();
    faceCentre_ = vector::zero;
    faceArea_ = vector::zero;
    subCellCentre_ = vector::zero;
    subCellVolume_ = -10;
    VOF_ = -10;
    cellStatus_ = -1;
    freeSurf_ = triSurface();
}

// ************************************************************************* //
