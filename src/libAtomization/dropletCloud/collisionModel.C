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

#include "dropletCloud.H"
#include "collisionModel.H"


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


bool Foam::collisionModel::collideDroplets
(
    const scalar dt,
    droplet& p1,
    droplet& p2,
    scalar& m1,
    scalar& m2
)
{
    bool coalescence = false;

    const vector pos1(p1.position());
    const vector pos2(p2.position());

    const vector& U1 = p1.U();
    const vector& U2 = p2.U();

    vector URel(U1 - U2);

    vector d(pos2 - pos1);
    scalar magd = mag(d);

    const scalar d1 = p1.d();
    const scalar d2 = p2.d();

    scalar sumD = d1 + d2;

// collision occurs
    if (0.5*sumD > magd)
    {
        coalescence = this->collideSorted(dt, p1, p2, m1, m2);
    }


    return coalescence;
}


bool Foam::collisionModel::collideSorted
(
    const scalar dt,
    droplet& p1,
    droplet& p2,
    scalar& m1,
    scalar& m2
)
{
    const scalar d1 = p1.d();
    const scalar d2 = p2.d();

    const vector& U1 = p1.U();
    const vector& U2 = p2.U();

    vector URel = U1 - U2;
    scalar magURel = mag(URel);

    scalar mTot = m1 + m2;

    p1.U() = (m1*U1 + m2*U2)/(m1+m2);
        
    m1 += m2;
    m2 = -1;

    return true;

}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


Foam::collisionModel::collisionModel
(
    const fvMesh& mesh
)
:
    mesh_(mesh),
    dict_
    (
        IOobject
        (
            "cloudProperties",
            mesh_.time().constant(),
            mesh_,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    ),
    active_(readBool(dict_.lookup("collision"))),
    ranGen_(clock::getTime() + pid()),
    cTime_(dict_.lookupOrDefault("cTime", 1.0)),
    cSpace_(dict_.lookupOrDefault("cSpace", 0.3)),
    rhop_(0.0),
    sigma_(0.0)
{
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void Foam::collisionModel::update
(
    dropletCloud& cloud,
    const scalar dt
)
{
    if (!active_)
    {
        return;
    }

    // Store cloud droplet properties
    rhop_ = cloud.rhop();
    sigma_ = cloud.sigma();

    // Create pointer to cloud
    dropletCloud *ptrCloud = const_cast<dropletCloud*>(&cloud);    

    // Create the occupancy list for the cells
    labelList occupancy(mesh_.nCells(), 0);
    forAllIter(Cloud<droplet>, *ptrCloud,  iter)
    {
        occupancy[iter().cell()]++;
    }

    // Initialize the sizes of the lists of parcels in each cell
    CompactListList<droplet*> pInCell(occupancy);

    // Reset the occupancy to use as a counter
    occupancy = 0;

    // Set the parcel pointer lists for each cell
    forAllIter(Cloud<droplet>, *ptrCloud, iter)
    {
        pInCell(iter().cell(), occupancy[iter().cell()]++) = &iter();
	droplet& p1 = iter();
	
	if (p1.d() > ROOTVSMALL)
	{
	forAllIter(Cloud<droplet>, *ptrCloud, iter2)
	{
	    droplet& p2 = iter2();
	    
	  if(p2.d() > ROOTVSMALL && p2 != p1)
	  {
	    if( 0.5*(p2.d() + p1.d()) > mag(p1.position() - p2.position()) )
        {
		  scalar m1 =  (4.0/3.0)*constant::mathematical::pi*pow(p1.d()/2.0, 3.0);
		  scalar m2 = (4.0/3.0)*constant::mathematical::pi*pow(p2.d()/2.0, 3.0);
		  
		  p1.U() = (m1*p1.U() + m2*p2.U())/(m1+m2);
		  m1 = m1 + m2;

		  p1.position() = 0.5*(p1.position() + p2.position());
		  p1.d() = cbrt(6.0*m1/(constant::mathematical::pi));
		  p2.d() = -1;
        }
			
	  }
	  }
	}
    }

    // Delete all droplets with negative diameter
    forAllIter(Cloud<droplet>, *ptrCloud, iter)
    {
        droplet& p = iter();
        
        if (p.d() < VSMALL)
        {
            cloud.deleteParticle(p);
        }
    }
}

// ************************************************************************* //
