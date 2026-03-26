/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
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

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inline Foam::droplet::trackingData::trackingData
(
    const dropletCloud& spc,
    const interpolationCellPoint<scalar>& rhoInterp,
    const interpolationCellPoint<vector>& UInterp,
    const interpolationCellPoint<scalar>& nuInterp,
    const vector& g
)
:
    particle::trackingData(spc),
    rhoInterp_(rhoInterp),
    UInterp_(UInterp),
    nuInterp_(nuInterp),
    g_(g)
{}


inline Foam::droplet::droplet
(
    const polyMesh& mesh,
    const vector& position,
    const label celli,
    const scalar d,
    const vector& U
)
:
    particle(mesh, position, celli),
    d_(d),
    U_(U),
    nParticle_(1),
    y_(0.0),
    yDot_(0.0)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline const Foam::interpolationCellPoint<Foam::scalar>&
Foam::droplet::trackingData::rhoInterp() const
{
    return rhoInterp_;
}


inline const Foam::interpolationCellPoint<Foam::vector>&
Foam::droplet::trackingData::UInterp() const
{
    return UInterp_;
}


inline const Foam::interpolationCellPoint<Foam::scalar>&
Foam::droplet::trackingData::nuInterp() const
{
    return nuInterp_;
}

inline const Foam::vector& Foam::droplet::trackingData::g() const
{
    return g_;
}


inline Foam::scalar& Foam::droplet::d()
{
    return d_;
}


inline Foam::vector& Foam::droplet::U()
{
    return U_;
}


inline Foam::scalar& Foam::droplet::nParticle()
{
    return nParticle_;
}


inline Foam::scalar& Foam::droplet::y()
{
    return y_;
}


inline Foam::scalar& Foam::droplet::yDot()
{
    return yDot_;
}


// ************************************************************************* //
