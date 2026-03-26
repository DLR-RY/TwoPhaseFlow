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
#include "fvMesh.H"
#include "volFields.H"
#include "interpolationCellPoint.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


void Foam::dropletCloud::createDropletDataFiles()
{
    // resize file pointer
    cloudDataFilesPtr_.resize(faceZoneIDs_.size());
    
    forAll (cloudDataFilesPtr_, i)
    {
        if (!cloudDataFilesPtr_[i] && Pstream::master() && faceZoneIDs_.size())
        {
            word startTimeName = this->db().time().timeName(this->db().time().startTime().value());

            fileName cloudDataFileDir = "postProcessing/dropletCloud/" + startTimeName + "/";

            fileName dataFile = faceZoneNames_[i] + "_cloudData.dat";
        
            mkDir(cloudDataFileDir);

            cloudDataFilesPtr_[i].reset
            (
                new OFstream
                (
                    cloudDataFileDir/dataFile
                )
            ); 

            cloudDataFilesPtr_[i]() << "// time" << tab << tab
                                    << "diameter" << tab << tab
                                    << "position" << tab << tab
                                    << "nParticle" << endl;
	}
    }
}


void Foam::dropletCloud::info()
{
    label cloudSize = (*this).size();
    reduce( cloudSize, sumOp<label>() );

    Info<< endl << "Droplet cloud: " << endl
        << "    Number of droplets : " << cloudSize << endl;
    if (cloudSize > 0)
    {
        Info<< "    Max diameter [mu]  : " << 1e6*Dmax() << endl
            << "    D10 diameter [mu]  : " << 1e6*Dij(1,0) << endl
            << "    D32 diameter [mu]  : " << 1e6*Dij(3,2) << endl;
    }
    Info<< endl;
}


Foam::scalar Foam::dropletCloud::Dij
(
    const label i,
    const label j
)
{
    scalar si = 0.0;
    scalar sj = 0.0;
    forAllIter(Cloud<droplet>, *this, iter)
    {
        droplet& p = iter();
        si += p.nParticle()*pow(p.d(), i);
        sj += p.nParticle()*pow(p.d(), j);
    }

    reduce(si, sumOp<scalar>());
    reduce(sj, sumOp<scalar>());

    return si/sj;
}


Foam::scalar Foam::dropletCloud::Dmax()
{
    scalar d = -GREAT;
    forAllIter(Cloud<droplet>, *this, iter)
    {
        droplet& p = iter();
        d = max(d, p.d());
    }

    reduce(d, maxOp<scalar>());

    return max(0.0, d);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


Foam::dropletCloud::dropletCloud
(
    const fvMesh& mesh,
    const dimensionedVector& g,
    const word& cloudName,
    bool readFields
)
:
    Cloud<droplet>(mesh, cloudName, false),
    mesh_(mesh),
    g_(g),
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
    phaseName_(word(dict_.lookup("phaseName"))),
    faceZoneIDs_(),
    faceZoneNames_(),
    cloudDataFilesPtr_(),
    dataDiameter_(),
    dataPosition_(),
    source_(mesh_.nCells(), Zero),
    momentumSource_
    (
        IOobject
        (
            "momentumSource",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedVector("zero", dimensionSet(1,-2,-2,0,0,0,0),vector::zero)
    ),
    collision_(mesh_),
    breakup_(mesh_)
{
    if (readFields)
    {
        droplet::readFields(*this);
    }

    // Collect face zones for postprocessing
    wordList faceZoneNames(dict_.lookup("faceZones"));
    const faceZoneMesh& fzm = mesh_.faceZones();
    DynamicList<label> zoneIDs;
    DynamicList<word> zoneNames;

    Info<< nl << type() << " faceZones" << endl;
    forAll(faceZoneNames, i)
    {
        const word& zoneName = faceZoneNames[i];
        label zoneI = fzm.findZoneID(zoneName);

        if (zoneI != -1)
        {
            zoneIDs.append(zoneI);
            zoneNames.append(zoneName);

	    const faceZone& fz = fzm[zoneI];
            label nFaces = returnReduce(fz.size(), sumOp<label>());
            Info<< "    " << zoneName << " faces: " << nFaces << nl;
	}
    }

    faceZoneNames_.transfer(zoneNames);
    faceZoneIDs_.transfer(zoneIDs);
    if(!faceZoneIDs_.size())
    {
        Info<< "    none" << endl;
    }
    Info << endl;

    // Get phase properties
    const dictionary& transportProperties = db().lookupObject<IOdictionary>
    (
        "transportProperties"
    );

    // surface tension
    rhop_ = readScalar(transportProperties.subDict(phaseName_).lookup("rho"));
    
    // dynamic viscosity
    mup_ = readScalar(transportProperties.subDict(phaseName_).lookup("nu"))*rhop_;

    // surface tension
    sigma_ = readScalar(transportProperties.lookup("sigma"));

    // Create postProcessing data file
    createDropletDataFiles();
    dataDiameter_.setSize(Pstream::nProcs());
    dataPosition_.setSize(Pstream::nProcs());
    dataNParticle_.setSize(Pstream::nProcs());
    forAll(dataDiameter_, i)
    {
        dataDiameter_[i].setSize(faceZoneIDs_.size());
        dataNParticle_[i].setSize(faceZoneIDs_.size());
        dataPosition_[i].setSize(faceZoneIDs_.size());
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void Foam::dropletCloud::move()
{
    // Reset data collection from processors
    dataDiameter_.clear();
    dataPosition_.clear();
    dataNParticle_.clear();
    dataDiameter_.setSize(Pstream::nProcs());
    dataPosition_.setSize(Pstream::nProcs());
    dataNParticle_.setSize(Pstream::nProcs());
    forAll(dataDiameter_, i)
    {
        dataDiameter_[i].setSize(faceZoneIDs_.size());
        dataPosition_[i].setSize(faceZoneIDs_.size());
        dataNParticle_[i].setSize(faceZoneIDs_.size());
    }

    const volScalarField& rho = mesh_.lookupObject<const volScalarField>("rho");
    const volVectorField& U = mesh_.lookupObject<const volVectorField>("U");
    const volScalarField& nu = mesh_.lookupObject<const volScalarField>("nu");

    interpolationCellPoint<scalar> rhoInterp(rho);
    interpolationCellPoint<vector> UInterp(U);
    interpolationCellPoint<scalar> nuInterp(nu);
        
    // Reset source term
    source_.resize(mesh_.nCells());
    source_ = vector::zero;

    droplet::trackingData
        td(*this, rhoInterp, UInterp, nuInterp, g_.value());

    // Move droplets
    Cloud<droplet>::move(*this, td, mesh_.time().deltaTValue());

    // Droplet collision
    collision_.update(*this, mesh_.time().deltaTValue());

    // Secondary breakup
    //breakup_.update(*this, td, mesh_.time().deltaTValue());

    // Source term for momentum equation
    //momentumSource_.primitiveFieldRef() = source_ / (mesh_.time().deltaT().value() * mesh_.V());

    // Collect data from all processors for postprocessing and write to file
    Pstream::gatherList(dataDiameter_);
    Pstream::gatherList(dataPosition_);
    Pstream::gatherList(dataNParticle_);
    for (int i=0; i<Pstream::nProcs(); i++)
    {
        // loop over all face zones
        forAll(dataDiameter_[i], j)
        {
            // loop over all data points
	    forAll(dataDiameter_[i][j], k)
            {
                if(Pstream::master())
                {
                    cloudDataFilesPtr_[j]() << this->db().time().value() << tab << dataDiameter_[i][j][k] << tab << dataPosition_[i][j][k] << tab << dataNParticle_[i][j][k]<< endl;
                }
            }
        }
    }

    info();
}


void Foam::dropletCloud::inject
(
    vector position,
    const scalar diameter,
    const vector velocity
)
{
    label celli;
    label tetFacei;
    label tetPti;

    volVectorField cellCentres = mesh_.C();

    // Find the cell, tetFacei and tetPti for point position
    mesh_.findCellFacePt
    (
        position,
        celli,
        tetFacei,
        tetPti
    );

    label proci = -1;

    if (celli >= 0)
    {
        proci = Pstream::myProcNo();
    }

    reduce(proci, maxOp<label>());

    // Ensure that only one processor attempts to insert this Parcel
    if (proci != Pstream::myProcNo())
    {
        celli = -1;
        tetFacei = -1;
        tetPti = -1;
    }

    // Last chance - find nearest cell and try that one - the point is
    // probably on an edge
    if (proci == -1)
    {
        celli = mesh_.findNearestCell(position);

        if (celli >= 0)
        {
            position += 1e-3*(cellCentres[celli] - position);

            mesh_.findCellFacePt
            (
                position,
                celli,
                tetFacei,
                tetPti
            );

            if (celli > 0)
            {
                proci = Pstream::myProcNo();
            }
        }

        reduce(proci, maxOp<label>());

        if (proci != Pstream::myProcNo())
        {
            celli = -1;
            tetFacei = -1;
            tetPti = -1;
        }
    }

    if (celli > -1)
    {
        droplet* pPtr = new droplet(mesh_, position, celli, diameter, velocity);
        
        Cloud<droplet>::addParticle(pPtr);
    }
}



// ************************************************************************* //
