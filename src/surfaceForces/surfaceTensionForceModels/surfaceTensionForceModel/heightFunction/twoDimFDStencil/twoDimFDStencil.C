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

#include "twoDimFDStencil.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::label Foam::twoDimFDStencil::arrayToList(const label& i,const label& j)
{
    // assumes that i and j are between 0 and 2
    return (i + 3*j);
}


Foam::label Foam::twoDimFDStencil::arrayToList3D
(
    const label& i,
    const label& j,
    const label& k
)
{
    // assumes that i,j and k are between 0 and 2
    return (i + 3*j + 9*k);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
Foam::twoDimFDStencil::twoDimFDStencil()
:
    twoDim_(false),
    heights_(9,0.0),
    dir_(0),
    gblcelli_(-1),
    status()
{

}


Foam::twoDimFDStencil::twoDimFDStencil
(
    const bool twoDim,
    const label dir,
    const label gblcelli
)
:
    twoDim_(twoDim),
    heights_(9),
    dir_(dir),
    gblcelli_(gblcelli),
    status()
{
    if (twoDim_)
    {
        heights_.setSize(3);
    }
    heights_ = 0;
    if (twoDim_ && dir >= 2)
    {
        FatalErrorInFunction << "the problem is two dimensional,"
                             << "direction has to be between 0 or 1"
                             << "dir is " << dir
                             << abort(FatalError);
    }
}


Foam::twoDimFDStencil::twoDimFDStencil
(
    const bool twoDim,
    const scalarField& heights,
    const label dir,
    const label gblcelli
)
:
    twoDim_(twoDim),
    heights_(heights),
    dir_(dir),
    gblcelli_(gblcelli),
    status()
{
    if (twoDim_ && dir >= 2)
    {
        FatalErrorInFunction << "the problem is two dimensional,"
                             << "direction has to be between 0 or 1"
                             << "dir is " << dir
                             << abort(FatalError);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar
Foam::twoDimFDStencil::calcCurvature(const scalar DeltaX)
{
    // http://basilisk.fr/src/curvature.h
    // the height are multiplied by the deltaX
    if (heights_.size() == 9 ) // 3D
    {
        scalar fx =
            (heights_[arrayToList(2,1)] - heights_[arrayToList(0,1)])/2;

        scalar fy =
            (heights_[arrayToList(1,2)] - heights_[arrayToList(1,0)])/2;

        scalar fxx =
        (
            heights_[arrayToList(2,1)]
        - 2*heights_[arrayToList(1,1)]
          + heights_[arrayToList(0,1)]
        )/(DeltaX);

        scalar fyy =
        (
            heights_[arrayToList(1,2)]
        - 2*heights_[arrayToList(1,1)]
          + heights_[arrayToList(1,0)]
        )/(DeltaX);

        scalar fxy =
        (
            heights_[arrayToList(2,2)]
          - heights_[arrayToList(2,0)]
          - heights_[arrayToList(0,2)]
          + heights_[arrayToList(0,0)]
        )/(4*DeltaX);

        // from jibben et al 2018
        return
           -(
                fxx + fyy
                + fxx*sqr(fy) + fyy*sqr(fx)
                - 2*fxy*fx*fy
            )/pow(1+sqr(fx)+sqr(fy),1.5);
    }
    else if (heights_.size() == 3)
    {

        scalar fx = (heights_[2]-heights_[0])/2;
        scalar fxx = (heights_[2]-2*heights_[1]+heights_[0])/(DeltaX);
        return -(fxx/pow(1+sqr(fx),1.5));
    }
    else
    {
        FatalErrorInFunction
           << " size should be 3 for the 2D case or 9 for the 3D case"
            << "the size is " << heights_.size()
            << abort(FatalError);
        return GREAT; // for the compiler
    }
}


Foam::scalar
Foam::twoDimFDStencil::addColumnHeight(const scalarList& stencilValues)
{
    // orientation true positive direction
    // orientation false negative direction
    scalar avgColumnHeight = 0;

   if (stencilValues.size() == 27) // 3D
   {
        if (dir_ == 0)
        {
            label i = 0,j = 0;
            forAll(heights_,hI)
            {
                scalar val = stencilValues[arrayToList3D(1,i,j)];
                avgColumnHeight += val;
                heights_[hI] += val;
                ++i;
                if (i == 3)
                {
                    i = 0;
                    ++j;
                }
            }

            return avgColumnHeight/heights_.size();
        }
        else if (dir_ == 1)
        {
            label i = 0,j = 0;
            forAll(heights_,hI)
            {
                scalar val = stencilValues[arrayToList3D(i,1,j)];
                avgColumnHeight += val;
                heights_[hI] += val;
                ++i;
                if (i == 3)
                {
                    i = 0;
                    ++j;
                }
            }

            return avgColumnHeight/heights_.size();

        }
        else if (dir_ == 2)
        {
            label i = 0,j = 0;
            forAll(heights_,hI)
            {
                scalar val = stencilValues[arrayToList3D(i,j,1)];
                avgColumnHeight += val;
                heights_[hI] += val;
                ++i;
                if (i == 3)
                {
                    i = 0;
                    ++j;
                }
            }

            return avgColumnHeight/heights_.size();

        }
        else
        {
            FatalErrorInFunction
                << "dir has to be between zero and 2 "
                << "dir_ : " << dir_
                << abort(FatalError);
        }
   }
   else if (stencilValues.size() == 9) // 2D case
   {
       if (dir_ == 0)
        {
            label i = 0;
            forAll(heights_,hI) // should be three
            {
                scalar val = stencilValues[arrayToList(1,i)];
                avgColumnHeight += val;
                heights_[hI] += val;
                i++;
            }

            return avgColumnHeight/heights_.size();

       }
       else if (dir_ == 1)
       {
            label i = 0;
            forAll(heights_,hI)
            {
                scalar val = stencilValues[arrayToList(i,1)];
                avgColumnHeight += val;
                heights_[hI] += val;
                i++;
            }

            return avgColumnHeight/heights_.size();

       }
   }
   else
   {
        FatalErrorInFunction
            << "stencil size should be 9 for the 2D case "
            << "or 27 for the 3D case the size is " << stencilValues.size()
            << abort(FatalError);
        return -1;
   }

    FatalErrorInFunction
        << "stencil size should be 9 for the 2D case "
        << "or 27 for the 3D case the size is " << stencilValues.size()
        << abort(FatalError);

    return -1;
}


Foam::Vector<Foam::label>
Foam::twoDimFDStencil::nextCell(HFStencil::orientation orientation)
{
    Vector<label> posInStencil(1,1,1);
    if (twoDim_)
    {
        posInStencil.z() = 0;
    }

    if (orientation == HFStencil::orientation::pos)
    {
        posInStencil[dir_] = 2;
    }
    else
    {
        posInStencil[dir_] = 0;
    }

    return posInStencil;
}


Foam::Ostream& Foam::operator<<(Ostream& os, const twoDimFDStencil& FDS)
{
    if (os.format() == IOstream::ASCII)
    {
        os  << FDS.twoDim_ << token::SPACE
            << FDS.heights_ << token::SPACE
            << FDS.dir_  << token::SPACE
            << FDS.gblcelli_ << token::SPACE
            << FDS.status[0].avgColVal  << token::SPACE
            << FDS.status[0].iterI  << token::SPACE
            << FDS.status[0].gblIdx << token::SPACE
            << FDS.status[1].avgColVal  << token::SPACE
            << FDS.status[1].iterI  << token::SPACE
            << FDS.status[1].gblIdx;
    }
    else
    {
        os  << FDS.twoDim_
            << FDS.heights_
            << FDS.dir_
            << FDS.gblcelli_
            << FDS.status[0].avgColVal
            << FDS.status[0].iterI
            << FDS.status[0].gblIdx
            << FDS.status[1].avgColVal
            << FDS.status[1].iterI
            << FDS.status[1].gblIdx;
    }


    return os;
}


Foam::Istream& Foam::operator>>(Istream& is, twoDimFDStencil& FDS)
{
    is  >> FDS.twoDim_
        >> FDS.heights_
        >> FDS.dir_
        >> FDS.gblcelli_
        >> FDS.status[0].avgColVal
        >> FDS.status[0].iterI
        >> FDS.status[0].gblIdx
        >> FDS.status[1].avgColVal
        >> FDS.status[1].iterI
        >> FDS.status[1].gblIdx;

    return is;
}


bool Foam::twoDimFDStencil::operator==
(
    const twoDimFDStencil& FDS2
) const
{
    return  (twoDim_ ==  FDS2.twoDim_) &&
            (heights_ ==  FDS2.heights_) &&
            (dir_ ==  FDS2.dir_);
}
// ************************************************************************* //
