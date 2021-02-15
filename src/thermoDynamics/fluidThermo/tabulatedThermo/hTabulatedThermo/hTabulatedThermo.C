/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2017 OpenFOAM Foundation
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

#include "hTabulatedThermo.H"
#include "IOstreams.H"


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

template<class EquationOfState>
void Foam::hTabulatedThermo<EquationOfState>::check
(
    List<Tuple2<scalar, scalar>>& data
) const
{
    const label n = data.size();
    scalar prevValue = data.first().first();

    for (label i=1; i<n; ++i)
    {
        const scalar currValue = data[i].first();

        // avoid duplicate values (divide-by-zero error)
        if (currValue <= prevValue)
        {
            FatalErrorInFunction
                << "out-of-order value: "
                << currValue << " at index " << i << nl
                << exit(FatalError);
        }
        prevValue = currValue;
    }
}


template<class EquationOfState>
void Foam::hTabulatedThermo<EquationOfState>::readTable
(
    tableReader<scalar>& reader,
    List<Tuple2<scalar, scalar>>& data,
    fileName file
)
{
    // preserve the original (unexpanded) fileName to avoid absolute paths
    // appearing subsequently in the write() method
    fileName fName(file);

    fName.expand();

    // Read data from file
    reader(fName, data);

    if (data.empty())
    {
        FatalErrorInFunction
            << "table read from " << fName << " is empty" << nl
            << exit(FatalError);
    }

    // Check that the data are okay
    check(data);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class EquationOfState>
Foam::hTabulatedThermo<EquationOfState>::hTabulatedThermo
(
    const dictionary& dict
)
:
    EquationOfState(dict),
    Hf_(dict.subDict("thermodynamics").get<scalar>("Hf")),
    Sf_(dict.subDict("thermodynamics").get<scalar>("Sf")),

    // fileName_(dict.subDict("thermodynamics").get<fileName>("CpFile")),
    // reader_(tableReader<scalar>::New(dict.subDict("thermodynamics"))),
    CpTable_()
{
    scalar nPoints = dict.subDict("transport").get<scalar>("nPoints");
    scalar minTVal = dict.subDict("thermodynamics").get<scalar>("TMin");
    scalar maxTVal = dict.subDict("thermodynamics").get<scalar>("TMax");
    fileName file = dict.subDict("thermodynamics").get<fileName>("CpFile");
    // hCoeffs_ = CpCoeffs_.integral();
    // sCoeffs_ = CpCoeffs_.integralMinus1();

    // // Offset h poly so that it is relative to the enthalpy at Tstd
    // hCoeffs_[0] += Hf_ - hCoeffs_.value(Tstd);

    // // Offset s poly so that it is relative to the entropy at Tstd
    // sCoeffs_[0] += Sf_ - sCoeffs_.value(Tstd);
    autoPtr<tableReader<scalar>> reader =
        tableReader<scalar>::New(dict.subDict("thermodynamics"));

    List<Tuple2<scalar, scalar>> data;
    readTable(reader.ref(),data,file);


    scalar deltaT = maxTVal - minTVal;
    scalar sample = (deltaT/nPoints);

    if (!data.empty())
    {
        CpTable_.first().setSize(nPoints);
        CpTable_.second().setSize(nPoints);
    }
    scalarField Temp(data.size());
    scalarField cp(data.size());

    forAll(data,i)
    {
        Temp[i] = data[i].first();
        cp[i] = data[i].second();
    }

    for (label n = 0;n < nPoints;n++)
    {
        scalar T = minTVal + n*sample;
        CpTable_.first()[n] = T;
        CpTable_.second()[n] = interpolateXY(T,Temp,cp);
    }

}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class EquationOfState>
void Foam::hTabulatedThermo<EquationOfState>::write
(
    Ostream& os
) const
{
    EquationOfState::write(os);

    // Entries in dictionary format
    {
        os.beginBlock("thermodynamics");
        os.writeEntry("Hf", Hf_);
        os.writeEntry("Sf", Sf_);
        //os.writeEntry(coeffsName("Cp"), CpCoeffs_);
        os.endBlock();
    }
}


// * * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * //

template<class EquationOfState>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const hTabulatedThermo<EquationOfState>& pt
)
{
    pt.write(os);
    return os;
}


// ************************************************************************* //
