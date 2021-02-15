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

#include "tabulatedTransport.H"
#include "IOstreams.H"
#include "tableReader.H"

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

template<class Thermo>
void Foam::tabulatedTransport<Thermo>::check
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


template<class Thermo>
void Foam::tabulatedTransport<Thermo>::readTable
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

template<class Thermo>
Foam::tabulatedTransport<Thermo>::tabulatedTransport
(
    const dictionary& dict
)
:
    Thermo(dict),
    // interpolMuTable_(dict.subDict("transport").get<fileName>("MuFile")),
    // interpolKappaTable_(dict.subDict("transport").get<fileName>("KappaFile"))
    MuTable_(),
    KappaTable_()
{
    scalar nPoints = dict.subDict("transport").get<scalar>("nPoints");
    scalar minTVal = dict.subDict("transport").get<scalar>("TMin");
    scalar maxTVal = dict.subDict("transport").get<scalar>("TMax");
    fileName muFile = dict.subDict("transport").get<fileName>("MuFile");
    fileName kappaFile = dict.subDict("transport").get<fileName>("KappaFile");

    autoPtr<tableReader<scalar>> reader =
        tableReader<scalar>::New(dict.subDict("transport"));

    List<Tuple2<scalar, scalar>> dataMu;
    List<Tuple2<scalar, scalar>> dataKappa;
    readTable(reader.ref(),dataMu,muFile);
    readTable(reader.ref(),dataKappa,kappaFile);


    if (!dataMu.empty())
    {
        MuTable_.first().setSize(nPoints+1);
        MuTable_.second().setSize(nPoints+1);
    }
    scalarField TempMu(dataMu.size());
    scalarField mu(dataMu.size());

    forAll(dataMu,i)
    {
        TempMu[i] = dataMu[i].first();
        mu[i] = dataMu[i].second();
    }

    if (!dataKappa.empty())
    {
        KappaTable_.first().setSize(nPoints+1);
        KappaTable_.second().setSize(nPoints+1);
    }
    scalarField TempKappa(dataKappa.size());
    scalarField kappa(dataKappa.size());

    forAll(dataKappa,i)
    {
        TempKappa[i] = dataKappa[i].first();
        kappa[i] = dataKappa[i].second();
    }

    scalar deltaT = maxTVal - minTVal;
    scalar sample = (deltaT/nPoints);

    for (label n = 0;n <= nPoints;n++)
    {
        scalar T = minTVal + n*sample;
        MuTable_.first()[n] = T;
        MuTable_.second()[n] = interpolateXY(T,TempMu,mu);
        KappaTable_.first()[n] = T;
        KappaTable_.second()[n] = interpolateXY(T,TempKappa,kappa);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Thermo>
void Foam::tabulatedTransport<Thermo>::write(Ostream& os) const
{
    os.beginBlock(this->name());

    Thermo::write(os);

    // Entries in dictionary format
    {
        os.beginBlock("transport");
        //os.writeEntry(coeffsName("mu"), muCoeffs_);
        //os.writeEntry(coeffsName("kappa"), kappaCoeffs_);
        os.endBlock();
    }

    os.endBlock();
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

template<class Thermo>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const tabulatedTransport<Thermo>& pt
)
{
    pt.write(os);
    return os;
}


// ************************************************************************* //
