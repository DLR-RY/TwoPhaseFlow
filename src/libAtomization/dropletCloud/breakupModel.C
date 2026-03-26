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
#include "breakupModel.H"


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


void Foam::breakupModel::breakup
(
    const scalar dt,
    const vector g,
    scalar& d,
    scalar& nParticle,
    scalar& y,
    scalar& yDot,
    const scalar rhop,
    const scalar mup,
    const scalar sigma,
    const vector& Up,
    const scalar rhoc,
    const scalar muc,
    const scalar Urmag
)
{
    if(model_ == "ETAB")
    {
        scalar r = 0.5*d;
        scalar r2 = r*r;
        scalar r3 = r*r2;

        scalar semiMass = nParticle*pow3(d);

        // Inverse of characteristic viscous damping time
        scalar rtd = 0.5*TABCmu_*mup/(rhop*r2);

        // Oscillation frequency (squared)
        scalar omega2 = TABComega_*sigma/(rhop*r3) - rtd*rtd;

        if (omega2 > 0)
        {
            scalar omega = sqrt(omega2);
            scalar romega = 1.0/omega;

            scalar We = rhoc*sqr(Urmag)*r/sigma;
            scalar Wetmp = We/TABtwoWeCrit_;

            // Initial values for y and yDot
            scalar y0 = y - We;
            scalar yDot0 = yDot + y0*rtd;

            // Update distortion parameters
            scalar c = cos(omega*dt);
            scalar s = sin(omega*dt);
            scalar e = exp(-rtd*dt);

            // Update distortion
            y = We + e*(y0*c + (yDot0/omega)*s);
            yDot = (We - y)*rtd + e*(yDot0*c - omega*y0*s);

            scalar y1 = y - Wetmp;
            scalar y2 = yDot*romega;

            scalar a = sqrt(y1*y1 + y2*y2);

            // scotty we may have break-up
            if (a + Wetmp > 1.0)
            {
                scalar phic = y1/a;

                // constrain phic within -1 to 1
                phic = max(min(phic, 1), -1);

                scalar phit = acos(phic);
                scalar phi = phit;
                scalar quad = -y2/a;
                if (quad < 0)
                {
                    phi = constant::mathematical::twoPi - phit;
                }

                scalar tb = 0;

                if (mag(y) < 1.0)
                {
                    scalar theta = acos((1.0 - Wetmp)/a);

                    if (theta < phi)
                    {
                        if (constant::mathematical::twoPi - theta >= phi)
                        {
                            theta = -theta;
                        }
                        theta += constant::mathematical::twoPi;
                    }
                    tb = (theta - phi)*romega;

                    // breakup occurs
                    if (dt > tb)
                    {
                        y = 1.0;
                        yDot = -a*omega*sin(omega*tb + phi);
                    }
                }
     
                // update droplet size
                if (dt > tb)
                {
                    scalar sqrtWe = AWe_*pow4(We) + 1.0;
                    scalar Kbr = k1_*omega*sqrtWe;

                    if (We > WeTransition_)
                    {
                        sqrtWe = sqrt(We);
                        Kbr =k2_*omega*sqrtWe;
                    }

                    scalar rWetmp = 1.0/Wetmp;
                    scalar cosdtbu = max(-1.0, min(1.0, 1.0 - rWetmp));
                    scalar dtbu = romega*acos(cosdtbu);
                    scalar decay = exp(-Kbr*dtbu);

                    scalar rNew = decay*r;
                    if (rNew < r)
                    {
                        d = 2.0*rNew;
                        y = 0.0;
                        yDot = 0.0;
                    }
                }
            }
        }
        else
        {
            // reset droplet distortion parameters
            y = 0;
            yDot = 0;
        }

        // update the nParticle count to conserve mass
        nParticle = semiMass/pow3(d);
    }
    else if (model_ == "ReitzDiwakar")
    {
        scalar d1 = d;
        scalar nuc = muc/rhoc;
        scalar We = 0.5*rhoc*sqr(Urmag)*d/sigma;
        scalar Re = Urmag*d/nuc;

        if (We > Cbag_)
        {
            if (We > Cstrip_*sqrt(Re))
            {
                scalar dStrip = sqr(2.0*Cstrip_*sigma)/(rhoc*pow3(Urmag)*muc);
                scalar tauStrip = Cs_*d*sqrt(rhop/rhoc)/Urmag;
                scalar fraction = dt/tauStrip;

                // new droplet diameter, implicit calculation
                d = (fraction*dStrip + d)/(1.0 + fraction);
            }
            else
            {
                scalar dBag = 2.0*Cbag_*sigma/(rhoc*sqr(Urmag));
                scalar tauBag = Cb_*d*sqrt(rhop*d/sigma);
                scalar fraction = dt/tauBag;

                // new droplet diameter, implicit calculation
                d = (fraction*dBag + d)/(1.0 + fraction);
            }

            // preserve the total mass/volume by updating the number of
            // particles in parcels due to breakup
            nParticle *= pow3(d1/d);
        }
    }
    else if (model_ == "PilchErdman")
    {
        // Weber number - eq (1)
        scalar We = rhoc*sqr(Urmag)*d/sigma;

        // Ohnesorge number - eq (2)
        scalar Oh = mup/sqrt(rhop*d*sigma);

        // Critical Weber number - eq (5)
        scalar Wec = 12.0*(1.0 + 1.077*pow(Oh, 1.6));

        if (We > Wec)
        {
            // We > 2670, wave crest stripping - eq (12)
            scalar taubBar = 5.5;

            if (We < 2670)
            {
                if (We > 351)
                {
                    // sheet stripping - eq (11)
                    taubBar = 0.766*pow(We - 12.0, 0.25);
                }
                else if (We > 45)
                {
                    // bag-and-stamen breakup  - eq (10)
                    taubBar = 14.1*pow(We - 12.0, 0.25);
                }
                else if (We > 18)
                {
                    // bag breakup - eq (9)
                    taubBar = 2.45*pow(We - 12.0, 0.25);
                }
                else if (We > 12)
                {
                    // vibrational breakup - eq (8)
                    taubBar = 6.0*pow(We - 12.0, -0.25);
                }
                else
                {
                    // no break-up
                    taubBar = GREAT;
                }
            }

            scalar rho12 = sqrt(rhoc/rhop);

            // velocity of fragmenting drop - eq (20)
            scalar Vd = Urmag*rho12*(B1_*taubBar + B2_*sqr(taubBar));

            // maximum stable diameter - eq (33)
            scalar Vd1 = sqr(1.0 - Vd/Urmag);
            Vd1 = max(Vd1, SMALL);
            scalar dStable = Wec*sigma/(Vd1*rhoc*sqr(Urmag));

            if (d < dStable)
            {
                // droplet diameter already stable = no break-up
                // - do not update d and nParticle
                return;
            }
            else
            {
                scalar semiMass = nParticle*pow3(d);

                // invert eq (3) to create a dimensional break-up time
                scalar taub = taubBar*d/(Urmag*rho12);

                // update droplet diameter according to the rate eq (implicitly)
                scalar frac = dt/taub;
                d = (d + frac*dStable)/(1.0 + frac);

                // correct the number of particles to conserve mass
                nParticle = semiMass/pow3(d);
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


Foam::breakupModel::breakupModel
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
    model_(word(dict_.lookup("breakupModel"))),
    k1_(0.2),
    k2_(0.2),
    WeTransition_(100.0),
    AWe_(0.0),
    TABComega_(8),
    TABCmu_(5),
    TABtwoWeCrit_(12),
    Cbag_(6.0),
    Cb_(0.785),
    Cstrip_(0.5),
    Cs_(10.0),
    B1_(0.375),
    B2_(0.2274)
{
    if
    (
        model_ != "PilchErdman"
     && model_ != "ETAB"
     && model_ != "ReitzDiwakar"
    )
    {
        Info << endl << "    Breakup model not found - disabled!" << endl;
        model_ = "none";
    }

    // ETAB constants
    scalar k21 = k2_/k1_;
    AWe_ = (k21*sqrt(WeTransition_) - 1.0)/pow4(WeTransition_);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void Foam::breakupModel::update
(
    dropletCloud& cloud,
    droplet::trackingData& td,
    const scalar dt
)
{
    if (model_ == "none")
    {
        return;
    }

    // Constant droplet properties
    const scalar rhop = cloud.rhop();
    const scalar mup = cloud.mup();
    const scalar sigma = cloud.sigma();

    // Create pointer to cloud
    dropletCloud *ptrCloud = const_cast<dropletCloud*>(&cloud);    

    // Loop over all droplets
    forAllIter(Cloud<droplet>, *ptrCloud,  iter)
    {
        // Cell values of carrier phase
        const tetIndices tetIs = iter().currentTetIndices();
        const scalar rhoc = td.rhoInterp().interpolate(iter().coordinates(), tetIs);
        const vector Uc = td.UInterp().interpolate(iter().coordinates(), tetIs);
        const scalar muc = td.nuInterp().interpolate(iter().coordinates(), tetIs) * rhoc;
	
        vector Urel = iter().U() - Uc; 
        const scalar Urmag = mag(Urel);

        const vector g = td.g();

        // Perform for breakup
        breakup
        (
            dt,
            g,
            iter().d(),
            iter().nParticle(),
            iter().y(),
            iter().yDot(),
            rhop,
            mup,
            sigma,
            iter().U(),
            rhoc,
            muc,
            Urmag
        );
    }
}

// ************************************************************************* //
