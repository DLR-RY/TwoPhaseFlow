/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2017 DHI
    Copyright (C) 2018-2019 Johan Roenby
    Copyright (C) 2019-2020 DLR
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

\*---------------------------------------------------------------------------*/

#include "cutFaceAdvect.H"
#include "OFstream.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::cutFaceAdvect::cutFaceAdvect
(
    const fvMesh& mesh,
    const volScalarField& alpha1
)
:
    cutFace(mesh),
    mesh_(mesh),
    alpha1_(alpha1),
    subFaceCentre_(Zero),
    subFaceArea_(Zero),
    subFacePoints_(10),
    surfacePoints_(4),
    pointStatus_(10),
    weight_(10),
    pTimes_(10),
    faceStatus_(-1)
{
    clearStorage();
}


// * * * * * * * * * * * Public Member Functions  * * * * * * * * * * * * * //

Foam::label Foam::cutFaceAdvect::calcSubFace
(
    const label faceI,
    const vector& normal,
    const vector& base
)
{
    clearStorage();

    const face& f = mesh_.faces()[faceI];

    label inLiquid = 0;
    label firstFullySubmergedPoint = -1;

    // Loop face and calculate pointStatus
    forAll(f, i)
    {
        scalar value = (mesh_.points()[f[i]] - base) & normal;
        if (mag(value) < 10 * SMALL)
        {
            value = 0;
        }
        pointStatus_.append(value);
        if (pointStatus_[i] > 10 * SMALL)
        {
            inLiquid++;
            if (firstFullySubmergedPoint == -1)
            {
                firstFullySubmergedPoint = i;
            }
        }
    }

    if (inLiquid == f.size()) // fluid face
    {
        faceStatus_ = -1;
        subFaceCentre_ = mesh_.faceCentres()[faceI];
        subFaceArea_ = mesh_.faceAreas()[faceI];
        return faceStatus_;
    }
    else if (inLiquid == 0) // gas face
    {
        faceStatus_ = 1;
        subFaceCentre_ = Zero;
        subFaceArea_ = Zero;
        return faceStatus_;
    }

    cutFace::calcSubFace
    (
        faceI,
        pointStatus_,
        firstFullySubmergedPoint,
        subFacePoints_,
        surfacePoints_,
        faceStatus_,
        subFaceCentre_,
        subFaceArea_
    );

    return faceStatus_;
}


Foam::label Foam::cutFaceAdvect::calcSubFace
(
    const face& f,
    const pointField& points,
    const scalarField& val,
    const scalar cutValue
)
{
    clearStorage();

    label inLiquid = 0;
    label firstFullySubmergedPoint = -1;
    scalarList pointStatus(f.size());

    // Loop face and calculate pointStatus
    forAll(f, i)
    {
        pointStatus[i] = val[f[i]] - cutValue;
        if (mag(pointStatus[i]) < 10 * SMALL)
        {
            pointStatus[i] = 0;
        }
        if (pointStatus[i] > 10 * SMALL)
        {
            ++inLiquid;
            if (firstFullySubmergedPoint == -1)
            {
                firstFullySubmergedPoint = i;
            }
        }
    }

    if (inLiquid == f.size()) // fluid face
    {
        faceStatus_ = -1;
        subFaceCentre_ = f.centre(points);
        subFaceArea_ = f.areaNormal(points);
        return faceStatus_;
    }
    else if (inLiquid == 0) // gas face
    {
        faceStatus_ = 1;
        subFaceCentre_ = Zero;
        subFaceArea_ = Zero;
        return faceStatus_;
    }

    cutFace::calcSubFace
    (
        f,
        points,
        pointStatus,
        firstFullySubmergedPoint,
        subFacePoints_,
        surfacePoints_,
        faceStatus_,
        subFaceCentre_,
        subFaceArea_
    );

    return faceStatus_;
}


Foam::scalar Foam::cutFaceAdvect::timeIntegratedFaceFlux
(
    const label faceI,
    const vector& x0,
    const vector& n0, // has to has the length of 1
    const scalar Un0,
    const scalar dt,
    const scalar phi,
    const scalar magSf
)
{
    clearStorage();

/* Temporarily taken out
    // Treating rare cases where isoface normal is not calculated properly
    if (mag(n0) < 0.5)
    {
        scalar alphaf = 0.0;
        scalar waterInUpwindCell = 0.0;

        if (phi > 0 || !mesh_.isInternalFace(faceI))
        {
            const label upwindCell = mesh_.faceOwner()[faceI];
            alphaf = alpha1_[upwindCell];
            waterInUpwindCell = alphaf * mesh_.V()[upwindCell];
        }
        else
        {
            const label upwindCell = mesh_.faceNeighbour()[faceI];
            alphaf = alpha1_[upwindCell];
            waterInUpwindCell = alphaf * mesh_.V()[upwindCell];
        }

        return min(alphaf * phi * dt, waterInUpwindCell);
    }*/

    // Find sorted list of times where the isoFace will arrive at face points
    // given initial position x0 and velocity Un0*n0

    // Get points for this face
    const face& f = mesh_.faces()[faceI];
    const label nPoints = f.size();

    if (mag(Un0) > 1e-12) // Note: tolerances
    {
        // Here we estimate time of arrival to the face points from their normal
        // distance to the initial surface and the surface normal velocity

        for (const scalar fi : f)
        {
            scalar value = ((mesh_.points()[fi] - x0) & n0) / Un0;
            if (mag(value) < 10 * SMALL)
            {
                value = 0;
            }
            pTimes_.append(value);
        }

        scalar dVf = 0;

        // Check if pTimes changes direction more than twice when looping face
        label nShifts = 0;
        forAll(pTimes_, pi) // i have no clue what this is
        {
            const label oldEdgeSign =
                sign(pTimes_[(pi + 1) % nPoints] - pTimes_[pi]);
            const label newEdgeSign =
                sign(pTimes_[(pi + 2) % nPoints] - pTimes_[(pi + 1) % nPoints]);

            if (newEdgeSign != oldEdgeSign)
            {
                nShifts++;
            }
        }

        if (nShifts == 2 || nShifts == 0)
        {
            dVf = phi / magSf * timeIntegratedArea(faceI, dt, magSf, Un0);
        }
        else if (nShifts >  2) // triangle decompose the non planar face
        {
            const pointField fPts(f.points(mesh_.points()));
            pointField fPts_tri(3);
            scalarField pTimes_tri(3);
            fPts_tri[0] = mesh_.faceCentres()[faceI];
            pTimes_tri[0] = ((fPts_tri[0] - x0) & n0)/Un0;
            scalar area = 0;
            for (label pi = 0; pi < nPoints; ++pi)
            {
                fPts_tri[1] = fPts[pi];
                pTimes_tri[1] = pTimes_[pi];
                fPts_tri[2] = fPts[(pi + 1) % nPoints];
                pTimes_tri[2] = pTimes_[(pi + 1) % nPoints];
                const scalar magSf_tri =
                    mag
                    (
                        0.5
                       *(fPts_tri[2] - fPts_tri[0])
                       ^(fPts_tri[1] - fPts_tri[0])
                    );
                area += magSf_tri;
                const scalar phi_tri = phi*magSf_tri/magSf;
                dVf +=
                    phi_tri/magSf_tri
                   *timeIntegratedArea
                    (
                        fPts_tri,
                        pTimes_tri,
                        dt,
                        magSf_tri,
                        Un0
                    );
            }
        }
        else
        {
            if (debug)
            {
                WarningInFunction
                    << "Warning: nShifts = " << nShifts << " on face "
                    << faceI << " with pTimes = " << pTimes_
                    << " owned by cell " << mesh_.faceOwner()[faceI]
                    << endl;
            }
        }

        return dVf;
    }
    else
    {
        // Un0 is almost zero and isoFace is treated as stationary
        calcSubFace(faceI, -n0, x0);
        const scalar alphaf = mag(subFaceArea() / magSf);

        if (debug)
        {
            WarningInFunction
                << "Un0 is almost zero (" << Un0
                << ") - calculating dVf on face " << faceI
                << " using subFaceFraction giving alphaf = " << alphaf
                << endl;
        }

        return phi * dt * alphaf;
    }
}


Foam::scalar Foam::cutFaceAdvect::timeIntegratedFaceFlux
(
    const label faceI,
    const scalarField& times,
    const scalar Un0,
    const scalar dt,
    const scalar phi,
    const scalar magSf
)
{
    clearStorage();

    label nPoints = times.size();

    {
        // Here we estimate time of arrival to the face points from their normal
        // distance to the initial surface and the surface normal velocity

        pTimes_.append(times);

        scalar dVf = 0;

        // Check if pTimes changes direction more than twice when looping face
        label nShifts = 0;
        forAll(pTimes_, pi) // i have no clue what this is
        {
            const label oldEdgeSign =
                sign(pTimes_[(pi + 1) % nPoints] - pTimes_[pi]);
            const label newEdgeSign =
                sign(pTimes_[(pi + 2) % nPoints] - pTimes_[(pi + 1) % nPoints]);

            if (newEdgeSign != oldEdgeSign)
            {
                ++nShifts;
            }
        }

        if (nShifts == 2)
        {
            dVf = phi/magSf*timeIntegratedArea(faceI, dt, magSf, Un0);
        }
        // not possible to decompose face
        return dVf;
    }
}


Foam::scalar Foam::cutFaceAdvect::timeIntegratedArea
(
    const pointField& fPts,
    const scalarField& pTimes,
    const scalar dt,
    const scalar magSf,
    const scalar Un0
)
{
    // Initialise time integrated area returned by this function
    scalar tIntArea = 0.0;

    // Finding ordering of vertex points
    labelList order(fPts.size());
    sortedOrder(pTimes, order);
    const scalar firstTime = pTimes[order.first()];
    const scalar lastTime = pTimes[order.last()];

    // Dealing with case where face is not cut by surface during time interval
    // [0,dt] because face was already passed by surface
    if (lastTime <= 0)
    {
        // If all face cuttings were in the past and cell is filling up (Un0>0)
        // then face must be full during whole time interval
        tIntArea = magSf * dt * pos0(Un0);
        return tIntArea;
    }

    // Dealing with case where face is not cut by surface during time interval
    // [0, dt] because dt is too small for surface to reach closest face point
    if (firstTime >= dt)
    {
        // If all cuttings are in the future but non of them within [0,dt] then
        // if cell is filling up (Un0 > 0) face must be empty during whole time
        // interval
        tIntArea = magSf * dt * (1 - pos0(Un0));
        return tIntArea;
    }

    // If we reach this point in the code some part of the face will be swept
    // during [tSmall, dt-tSmall]. However, it may be the case that there are no
    // vertex times within the interval. This will happen sometimes for small
    // time steps where both the initial and the final face-interface
    // intersection line (FIIL) will be along the same two edges.

    // Face-interface intersection line (FIIL) to be swept across face
    DynamicList<point> FIIL(3);
    // Submerged area at beginning of each sub time interval time
    scalar initialArea = 0.0;
    // Running time keeper variable for the integration process
    scalar time = 0.0;

    // Special treatment of first sub time interval
    if (firstTime > 0)
    {
        // If firstTime > 0 the face is uncut in the time interval
        // [0, firstTime] and hence fully submerged in fluid A or B.
        // If Un0 > 0 cell is filling up and it must initially be empty.
        // If Un0 < 0 cell must initially be full(y immersed in fluid A).
        time = firstTime;
        initialArea = magSf * (1.0 - pos0(Un0));
        tIntArea = initialArea * time;
        cutPoints(fPts, pTimes, time, FIIL);
    }
    else
    {
        // If firstTime <= 0 then face is initially cut and we must
        // calculate the initial submerged area and FIIL:
        time = 0.0;
        // Note: calcSubFace assumes well-defined 2-point FIIL!!!!
        // calcSubFace(fPts, -sign(Un0)*pTimes, time);
        // calcSubFace(fPts, -sign(Un0)*pTimes, time)
        calcSubFace(face(identity(pTimes.size())), fPts, pTimes, time);
        initialArea = mag(subFaceArea());
        cutPoints(fPts, pTimes, time, FIIL);
    }

    // Making sorted array of all vertex times that are between max(0,firstTime)
    // and dt and further than tSmall from the previous time.
    DynamicList<scalar> sortedTimes(pTimes.size());
    {
        scalar prevTime = time;
        const scalar tSmall = max(1e-6*dt, 10*SMALL);

        for (const scalar timeI : order)
        {
            if (timeI > prevTime + tSmall && timeI <= dt)
            {
                sortedTimes.append(timeI);
                prevTime = timeI;
            }
        }
    }

    // Sweeping all quadrilaterals corresponding to the intervals defined above
    for (const scalar newTime : sortedTimes)
    {
        // New face-interface intersection line
        DynamicList<point> newFIIL(3);
        cutPoints(fPts, pTimes, newTime, newFIIL);

        // quadrilateral area coefficients
        scalar alpha = 0, beta = 0;
        quadAreaCoeffs(FIIL, newFIIL, alpha, beta);
        // Integration of area(t) = A*t^2+B*t from t = 0 to 1
        tIntArea += (newTime - time) *
                    (initialArea + sign(Un0) * (alpha/3.0 + 0.5*beta));
        // Adding quad area to submerged area
        initialArea += sign(Un0)*(alpha + beta);

        FIIL = newFIIL;
        time = newTime;
    }

    // Special treatment of last time interval
    if (lastTime > dt)
    {
        // FIIL will end up cutting the face at dt
        // New face-interface intersection line
        DynamicList<point> newFIIL(3);
        cutPoints(fPts, pTimes, dt, newFIIL);

        // quadrilateral area coefficients
        scalar alpha = 0, beta = 0;
        quadAreaCoeffs(FIIL, newFIIL, alpha, beta);
        // Integration of area(t) = A*t^2+B*t from t = 0 to 1
        tIntArea += (dt - time) *
                    (initialArea + sign(Un0)*(alpha / 3.0 + 0.5 * beta));
    }
    else
    {
        // FIIL will leave the face at lastTime and face will be fully in fluid
        // A or fluid B in the time interval from lastTime to dt.
        tIntArea += magSf*(dt - lastTime)*pos0(Un0);
    }

    return tIntArea;
}


void Foam::cutFaceAdvect::isoFacesToFile
(
    const DynamicList<List<point>>& faces,
    const word& filNam,
    const word& filDir
) const
{
    // Writing isofaces to vtk file for inspection in paraview

    fileName outputFile(filDir/(filNam + ".vtk"));

    mkDir(filDir);
    Info<< "Writing file: " << outputFile << endl;

    OFstream os(outputFile);
    os  << "# vtk DataFile Version 2.0" << nl
        << filNam << nl
        << "ASCII" << nl
        << "DATASET POLYDATA" << nl;

    label nPoints{0};
    for (const List<point>& f : faces)
    {
        nPoints += f.size();
    }

    os << "POINTS " << nPoints << " float" << nl;
    for (const List<point>& f : faces)
    {
        for (const point& p : f)
        {
            os << p.x() << ' ' << p.y() << ' ' << p.z() << nl;
        }
    }

    os  << "POLYGONS "
        << faces.size() << ' ' << (nPoints + faces.size()) << nl;

    label pointi = 0;
    for (const List<point>& f : faces)
    {
        label endp = f.size();
        os << endp;

        endp += pointi;

        while (pointi < endp)
        {
            os << ' ' << pointi;
            ++pointi;
        }
        os << nl;
    }
}


Foam::label Foam::cutFaceAdvect::calcSubFace
(
    const label faceI,
    const label sign,
    const scalar cutValue
)
{
    // clearStorage();
    const face& f = mesh_.faces()[faceI];
    label inLiquid = 0;
    label firstFullySubmergedPoint = -1;

    // Loop face and calculate pointStatus
    forAll(f, i)
    {
        scalar value = (sign * pTimes_[i] - cutValue);

        if (mag(value) < 10 * SMALL)
        {
            value = 0;
        }
        pointStatus_.append(value);
        if (pointStatus_[i] > 10 * SMALL)
        {
            inLiquid++;
            if (firstFullySubmergedPoint == -1)
            {
                firstFullySubmergedPoint = i;
            }
        }
    }

    if (inLiquid == f.size()) // fluid face
    {
        faceStatus_ = -1;
        subFaceCentre_ = mesh_.faceCentres()[faceI];
        subFaceArea_ = mesh_.faceAreas()[faceI];
        return faceStatus_;
    }
    else if (inLiquid == 0) // gas face
    {
        faceStatus_ = 1;
        subFaceCentre_ = Zero;
        subFaceArea_ = Zero;
        return faceStatus_;
    }

    cutFace::calcSubFace
    (
        faceI,
        pointStatus_,
        firstFullySubmergedPoint,
        subFacePoints_,
        surfacePoints_,
        faceStatus_,
        subFaceCentre_,
        subFaceArea_
    );

    return faceStatus_;
}


Foam::scalar Foam::cutFaceAdvect::timeIntegratedArea
(
    const label faceI,
    const scalar dt,
    const scalar magSf,
    const scalar Un0
)
{
    // Initialise time integrated area returned by this function
    scalar tIntArea = 0.0;

    // Finding ordering of vertex points
    labelList order(pTimes_.size());
    sortedOrder(pTimes_, order);
    const scalar firstTime = pTimes_[order.first()];
    const scalar lastTime = pTimes_[order.last()];

    // Dealing with case where face is not cut by surface during time interval
    // [0,dt] because face was already passed by surface
    if (lastTime <= 0)
    {
        // If all face cuttings were in the past and cell is filling up (Un0>0)
        // then face must be full during whole time interval
        tIntArea = magSf* dt * pos0(Un0);
        return tIntArea;
    }

    // Dealing with case where face is not cut by surface during time interval
    // [0, dt] because dt is too small for surface to reach closest face point
    if (firstTime >= dt)
    {
        // If all cuttings are in the future but non of them within [0,dt] then
        // if cell is filling up (Un0 > 0) face must be empty during whole time
        // interval
        tIntArea = magSf * dt * (1 - pos0(Un0));
        return tIntArea;
    }

    // If we reach this point in the code some part of the face will be swept
    // during [tSmall, dt-tSmall]. However, it may be the case that there are no
    // vertex times within the interval. This will happen sometimes for small
    // time steps where both the initial and the final face-interface
    // intersection line (FIIL) will be along the same two edges.

    // Face-interface intersection line (FIIL) to be swept across face
    DynamicList<point> FIIL(3);
    // Submerged area at beginning of each sub time interval time
    scalar initialArea = 0.0;
    // Running time keeper variable for the integration process
    scalar time = 0.0;

    // Special treatment of first sub time interval
    if (firstTime > 0)
    {
        // If firstTime > 0 the face is uncut in the time interval
        // [0, firstTime] and hence fully submerged in fluid A or B.
        // If Un0 > 0 cell is filling up and it must initially be empty.
        // If Un0 < 0 cell must initially be full(y immersed in fluid A).
        time = firstTime;
        initialArea = magSf * (1.0 - pos0(Un0));
        tIntArea = initialArea * time;
        cutPoints(faceI, time, FIIL);
    }
    else
    {
        // If firstTime <= 0 then face is initially cut and we must
        // calculate the initial submerged area and FIIL:
        time = 0.0;
        // Note: calcSubFace assumes well-defined 2-point FIIL!!!!
        calcSubFace(faceI, -sign(Un0), time);
        initialArea = mag(subFaceArea());
        cutPoints(faceI, time, FIIL);
    }

    // Making sorted array of all vertex times that are between max(0,firstTime)
    // and dt and further than tSmall from the previous time.
    DynamicList<scalar> sortedTimes(pTimes_.size());
    {
        scalar prevTime = time;
        const scalar tSmall = max(1e-6*dt, 10*SMALL);
        for (const label oI : order)
        {
            const scalar timeI = pTimes_[oI];
            if (timeI > prevTime + tSmall && timeI <= dt)
            {
                sortedTimes.append(timeI);
                prevTime = timeI;
            }
        }
    }

    // Sweeping all quadrilaterals corresponding to the intervals defined above
    for (const scalar newTime : sortedTimes)
    {
        // New face-interface intersection line
        DynamicList<point> newFIIL(3);
        cutPoints(faceI, newTime, newFIIL);

        // quadrilateral area coefficients
        scalar alpha = 0, beta = 0;

        quadAreaCoeffs(FIIL, newFIIL, alpha, beta);
        // Integration of area(t) = A*t^2+B*t from t = 0 to 1
        tIntArea +=
            (newTime - time)
          * (initialArea + sign(Un0)
          * (alpha / 3.0 + 0.5 * beta));
        // Adding quad area to submerged area
        initialArea += sign(Un0) * (alpha + beta);

        FIIL = newFIIL;
        time = newTime;
    }

    // Special treatment of last time interval
    if (lastTime > dt)
    {
        // FIIL will end up cutting the face at dt
        // New face-interface intersection line
        DynamicList<point> newFIIL(3);
        cutPoints(faceI, dt, newFIIL);

        // quadrilateral area coefficients
        scalar alpha = 0, beta = 0;
        quadAreaCoeffs(FIIL, newFIIL, alpha, beta);
        // Integration of area(t) = A*t^2+B*t from t = 0 to 1
        tIntArea +=
            (dt - time)
          * (initialArea + sign(Un0) * (alpha / 3.0 + 0.5 * beta));
    }
    else
    {
        // FIIL will leave the face at lastTime and face will be fully in fluid
        // A or fluid B in the time interval from lastTime to dt.
        tIntArea += magSf * (dt - lastTime) * pos0(Un0);
    }

    return tIntArea;
}


void Foam::cutFaceAdvect::quadAreaCoeffs
(
    const DynamicList<point>& pf0,
    const DynamicList<point>& pf1,
    scalar& alpha,
    scalar& beta
) const
{
    // Number of points in provided face-interface intersection lines
    const label np0 = pf0.size();
    const label np1 = pf1.size();

    // quad area coeffs such that area(t) = alpha*t^2 + beta*t.
    // With time interval normalised, we have full quadArea = alpha + beta
    // and time integrated quad area = alpha/3 + beta/2;
    alpha = 0.0;
    beta = 0.0;

    if (np0 && np1)
    {
        // Initialising quadrilateral vertices A, B, C and D
        vector A(pf0[0]);
        vector C(pf1[0]);
        vector B(pf0[0]);
        vector D(pf1[0]);

        if (np0 == 2)
        {
            B = pf0[1];
        }
        else if (np0 > 2)
        {
            WarningInFunction << "Vertex face was cut at pf0 = " << pf0 << endl;
        }

        if (np1 == 2)
        {
            D = pf1[1];
        }
        else if (np1 > 2)
        {
            WarningInFunction << "Vertex face was cut at pf1 = " << pf1 << endl;
        }

        // Swapping pf1 points if pf0 and pf1 point in same general direction
        // (because we want a quadrilateral ABCD where pf0 = AB and pf1 = CD)
        if (((B - A) & (D - C)) > 0)
        {
            vector tmp = D;
            D = C;
            C = tmp;
        }

        // Defining local coordinates (xhat, yhat) for area integration of swept
        // quadrilateral ABCD such that A = (0,0), B = (Bx,0), C = (Cx,Cy) and
        // D = (Dx,Dy) with Cy = 0 and Dy > 0.

        const scalar Bx = mag(B - A);

        vector xhat(Zero);
        if (Bx > 10 * SMALL)
        {
            // If |AB| > 0 ABCD we use AB to define xhat
            xhat = (B - A) / mag(B - A);
        }
        else if (mag(C - D) > 10 * SMALL)
        {
            // If |AB| ~ 0 ABCD is a triangle ACD and we use CD for xhat
            xhat = (C - D) / mag(C - D);
        }
        else
        {
            return;
        }

        // Defining vertical axis in local coordinates
        vector yhat = D - A;
        yhat -= ((yhat & xhat) * xhat);

        if (mag(yhat) > 10 * SMALL)
        {
            yhat /= mag(yhat);

            const scalar Cx = (C - A) & xhat;
            const scalar Cy = mag((C - A) & yhat);
            const scalar Dx = (D - A) & xhat;
            const scalar Dy = mag((D - A) & yhat);

            // area = ((Cx - Bx)*Dy - Dx*Cy)/6.0 + 0.25*Bx*(Dy + Cy);
            alpha = 0.5 * ((Cx - Bx) * Dy - Dx * Cy);
            beta = 0.5 * Bx * (Dy + Cy);
        }
    }
    else
    {
        WarningInFunction
            << "Vertex face was cut at " << pf0 << " and at "
            << pf1 << endl;
    }
}


void Foam::cutFaceAdvect::cutPoints
(
    const label faceI,
    const scalar f0,
    DynamicList<point>& cutPoints
)
{
    const face& f = mesh_.faces()[faceI];
    const label nPoints = f.size();
    scalar f1(pTimes_[0]);

    // Snapping vertex value to f0 if very close (needed for 2D cases)
    if (mag(f1 - f0) < 10 * SMALL)
    {
        f1 = f0;
    }

    forAll(f, pi)
    {
        label pi2 = (pi + 1) % nPoints;
        scalar f2 = pTimes_[pi2];

        // Snapping vertex value
        if (mag(f2 - f0) < 10 * SMALL)
        {
            f2 = f0;
        }

        if ((f1 < f0 && f2 > f0) || (f1 > f0 && f2 < f0))
        {
            const scalar s = (f0 - f1) / (f2 - f1);
            cutPoints.append
            (
                mesh_.points()[f[pi]]
              + s*(mesh_.points()[f[pi2]] - mesh_.points()[f[pi]])
            );
        }
        else if (f1 == f0)
        {
            cutPoints.append(mesh_.points()[f[pi]]);
        }
        f1 = f2;
    }

    if (cutPoints.size() > 2)
    {
        WarningInFunction
            << "cutPoints = " << cutPoints
            << " for pts = " << f.points(mesh_.points())
            << ", f - f0 = " << f - f0 << " and f0 = " << f0
            << endl;
    }
}


void Foam::cutFaceAdvect::cutPoints
(
    const pointField& pts,
    const scalarField& f,
    const scalar f0,
    DynamicList<point>& cutPoints
)
{
    const label nPoints = pts.size();
    scalar f1(f[0]);

    // Snapping vertex value to f0 if very close (needed for 2D cases)
    if (mag(f1 - f0) < 10 * SMALL)
    {
        f1 = f0;
    }

    forAll(pts, pi)
    {
        label pi2 = (pi + 1) % nPoints;
        scalar f2 = f[pi2];

        // Snapping vertex value
        if (mag(f2 - f0) < 10 * SMALL)
        {
            f2 = f0;
        }

        if ((f1 < f0 && f2 > f0) || (f1 > f0 && f2 < f0))
        {
            const scalar s = (f0 - f1) / (f2 - f1);
            cutPoints.append(pts[pi] + s * (pts[pi2] - pts[pi]));
        }
        else if (f1 == f0)
        {
            cutPoints.append(pts[pi]);
        }
        f1 = f2;
    }

    if (cutPoints.size() > 2)
    {
        WarningInFunction
            << "cutPoints = " << cutPoints << " for pts = " << pts
            << ", f - f0 = " << f - f0 << " and f0 = " << f0
            << endl;
    }
}


const Foam::point& Foam::cutFaceAdvect::subFaceCentre() const
{
    return subFaceCentre_;
}


const Foam::vector& Foam::cutFaceAdvect::subFaceArea() const
{
    return subFaceArea_;
}


const Foam::DynamicList<Foam::point>& Foam::cutFaceAdvect::subFacePoints() const
{
    return subFacePoints_;
}


const Foam::DynamicList<Foam::point>& Foam::cutFaceAdvect::surfacePoints() const
{
    return surfacePoints_;
}


void Foam::cutFaceAdvect::clearStorage()
{
    subFaceCentre_ = Zero;
    subFaceArea_ = Zero;
    subFacePoints_.clear();
    surfacePoints_.clear();
    pointStatus_.clear();
    pTimes_.clear();
    weight_.clear();
    faceStatus_ = -1;
}


// ************************************************************************* //
