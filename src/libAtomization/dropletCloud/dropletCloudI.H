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

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline const Foam::fvMesh& Foam::dropletCloud::mesh() const
{
    return mesh_;
}

inline Foam::scalar Foam::dropletCloud::rhop() const
{
    return rhop_;
}

inline Foam::scalar Foam::dropletCloud::mup() const
{
    return mup_;
}

inline Foam::scalar Foam::dropletCloud::sigma() const
{
    return sigma_;
}

inline Foam::vectorField& Foam::dropletCloud::source()
{
    return source_;
}

inline const Foam::volVectorField& Foam::dropletCloud::momentumSource()
{
    return momentumSource_;
}

inline const Foam::labelList& Foam::dropletCloud::faceZoneIDs()
{
    return faceZoneIDs_;
}

inline Foam::List<Foam::List<Foam::DynamicList<Foam::scalar> > >& Foam::dropletCloud::dataDiameter()
{
    return dataDiameter_;
}

inline Foam::List<Foam::List<Foam::DynamicList<Foam::scalar> > >& Foam::dropletCloud::dataNParticle()
{
    return dataNParticle_;
}

inline Foam::List<Foam::List<Foam::DynamicList<Foam::vector> > >& Foam::dropletCloud::dataPosition()
{
    return dataPosition_;
}

// ************************************************************************* //
