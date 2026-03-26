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

#include "droplet.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const std::size_t Foam::droplet::sizeofFields
(
    sizeof(droplet) - sizeof(particle)
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::droplet::droplet
(
    const polyMesh& mesh,
    Istream& is,
    bool readFields,
    bool newFormat
)
:
    particle(mesh, is, readFields, newFormat)
{
    if (readFields)
    {
        if (is.format() == IOstream::ASCII)
        {
            d_ = readScalar(is);
            is >> U_;
            nParticle_ = readScalar(is);
            y_ = readScalar(is);
            yDot_ = readScalar(is);
        }
        else
        {
            is.read(reinterpret_cast<char*>(&d_), sizeofFields);
        }
    }

    is.check(FUNCTION_NAME);
}


void Foam::droplet::readFields(Cloud<droplet>& c)
{
    bool valid = c.size();

    particle::readFields(c);

    IOField<scalar> d(c.fieldIOobject("d", IOobject::MUST_READ), valid);
    c.checkFieldIOobject(c, d);

    IOField<vector> U(c.fieldIOobject("U", IOobject::MUST_READ), valid);
    c.checkFieldIOobject(c, U);

    IOField<scalar> nParticle(c.fieldIOobject("nParticle", IOobject::MUST_READ), valid);
    c.checkFieldIOobject(c, nParticle);

    IOField<scalar> y(c.fieldIOobject("y", IOobject::MUST_READ), valid);
    c.checkFieldIOobject(c, y);

    IOField<scalar> yDot(c.fieldIOobject("yDot", IOobject::MUST_READ), valid);
    c.checkFieldIOobject(c, yDot);

    label i = 0;
    forAllIter(Cloud<droplet>, c, iter)
    {
        droplet& p = iter();

        p.d_ = d[i];
        p.U_ = U[i];
        p.nParticle_ = nParticle[i];
        p.y_ = y[i];
        p.yDot_ = yDot[i];
        i++;
    }
}


void Foam::droplet::writeFields(const Cloud<droplet>& c)
{
    particle::writeFields(c);

    label np = c.size();

    IOField<scalar> d(c.fieldIOobject("d", IOobject::NO_READ), np);
    IOField<vector> U(c.fieldIOobject("U", IOobject::NO_READ), np);
    IOField<scalar> nParticle(c.fieldIOobject("nParticle", IOobject::NO_READ), np);
    IOField<scalar> y(c.fieldIOobject("y", IOobject::NO_READ), np);
    IOField<scalar> yDot(c.fieldIOobject("yDot", IOobject::NO_READ), np);

    label i = 0;
    forAllConstIter(Cloud<droplet>, c, iter)
    {
        const droplet& p = iter();

        d[i] = p.d_;
        U[i] = p.U_;
        nParticle[i] = p.nParticle_;
        y[i] = p.y_;
        yDot[i] = p.yDot_;
        i++;
    }

    d.write(np > 0);
    U.write(np > 0);
    nParticle.write(np > 0);
    y.write(np > 0);
    yDot.write(np > 0);
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<(Ostream& os, const droplet& p)
{
    if (os.format() == IOstream::ASCII)
    {
        os  << static_cast<const particle&>(p)
            << token::SPACE << p.d_
            << token::SPACE << p.U_
            << token::SPACE << p.nParticle_;
    }
    else
    {
        os  << static_cast<const particle&>(p);
        os.write
        (
            reinterpret_cast<const char*>(&p.d_),
            droplet::sizeofFields
        );
    }

    os.check(FUNCTION_NAME);
    return os;
}


// ************************************************************************* //
