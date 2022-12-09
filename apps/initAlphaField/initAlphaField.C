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

Application
    initAlphaField

Description
    Uses cellCellIso to create a volume fraction field from either a cylinder,
    a sphere or a plane.

Author
    Henning Scheufler, DLR, all rights reserved.
    Johan Roenby, DHI, all rights reserved.


\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "searchableCylinder.H"
#include "searchableSurface.H"
#include "mathematicalConstants.H"
#include "SortableList.H"

#include "triSurface.H"
#include "triSurfaceTools.H"

#include "cutFace.H"
#include "cutCell.H"

#include "implicitFunction.H"
#include "cutCellIso.H"
#include "cutCellImpFunc.H"



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void isoFacesToFile
(
    const DynamicList< List<point> >& faces,
    const word filNam,
    const word filDir
)
{
    //Writing isofaces to vtk file for inspection in paraview

    mkDir(filDir);
    autoPtr<OFstream> vtkFilePtr;
    if (Pstream::parRun())
    {
        // Collect points from all the processors
        List<DynamicList<List<point>>> allProcFaces(Pstream::nProcs());
        allProcFaces[Pstream::myProcNo()] = faces;
        Pstream::gatherList(allProcFaces);

        if (Pstream::master())
        {
            Info << "Writing file: " << (filDir + "/" + filNam + ".vtk") << endl;
            vtkFilePtr.reset(new OFstream(filDir + "/" + filNam + ".vtk"));
            vtkFilePtr() << "# vtk DataFile Version 2.0" << endl;
            vtkFilePtr() << filNam << endl;
            vtkFilePtr() << "ASCII" << endl;
            vtkFilePtr() << "DATASET POLYDATA" << endl;

            face f;
            label nPoints(0);
            label fSize = 0;
            forAll(allProcFaces, proci)
            {
                const DynamicList<List<point>>& procFaces = allProcFaces[proci];


                forAll(procFaces,fi)
                {
                    nPoints += procFaces[fi].size();
                    fSize++ ;
                }

            }

            vtkFilePtr() << "POINTS " << nPoints << " float" << endl;

            forAll(allProcFaces, proci)
            {
                const DynamicList<List<point>>& procFaces = allProcFaces[proci];

                forAll(procFaces,fi)
                {
                    List<point> pf = procFaces[fi];
                    forAll(pf,pi)
                    {
                        point p = pf[pi];
                        vtkFilePtr() << p[0] << " " << p[1] << " " << p[2] << endl;
                    }
                }

            }

            vtkFilePtr() << "POLYGONS " << fSize << " " << nPoints + fSize << endl;

            label np = 0;

            forAll(allProcFaces, proci)
            {
                const DynamicList<List<point>>& procFaces = allProcFaces[proci];


                forAll(procFaces,fi)
                {
                    nPoints = procFaces[fi].size();
                    vtkFilePtr() << nPoints;
                    for (label pi = np; pi < np + nPoints; pi++ )
                    {
                        vtkFilePtr() << " " << pi;
                    }
                    vtkFilePtr() << "" << endl;
                    np += nPoints;
                }

            }
        }
    }
    else
    {
        Info << "Writing file: " << (filDir + "/" + filNam + ".vtk") << endl;
        vtkFilePtr.reset(new OFstream(filDir + "/" + filNam + ".vtk"));
        vtkFilePtr() << "# vtk DataFile Version 2.0" << endl;
        vtkFilePtr() << filNam << endl;
        vtkFilePtr() << "ASCII" << endl;
        vtkFilePtr() << "DATASET POLYDATA" << endl;
        label nPoints(0);
        forAll(faces,fi)
        {
            nPoints += faces[fi].size();
        }

        vtkFilePtr() << "POINTS " << nPoints << " float" << endl;
        forAll(faces,fi)
        {
            List<point> pf = faces[fi];
            forAll(pf,pi)
            {
                point p = pf[pi];
                vtkFilePtr() << p[0] << " " << p[1] << " " << p[2] << endl;
            }
        }
        vtkFilePtr() << "POLYGONS " << faces.size() << " " << nPoints + faces.size() << endl;

        label np = 0;
        forAll(faces,fi)
        {
            nPoints = faces[fi].size();
            vtkFilePtr() << nPoints;
            for (label pi = np; pi < np + nPoints; pi++ )
            {
                vtkFilePtr() << " " << pi;
            }
            vtkFilePtr() << "" << endl;
            np += nPoints;
        }
    }
}


int main(int argc, char *argv[])
{
    #include "addRegionOption.H"
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createNamedMesh.H"

    IOdictionary initAlphaFieldDict
    (
        IOobject
        (
            "initAlphaFieldDict",
            runTime.system(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    Info<< "Reading field alpha1\n" << endl;
    const word nameField (initAlphaFieldDict.lookup("field"));
    const bool invert = initAlphaFieldDict.lookupOrDefault<bool>("invert",false);
    const bool writeVTK = initAlphaFieldDict.lookupOrDefault<bool>("writeVTK",true);;

    volScalarField alpha1
    (
        IOobject
        (
            nameField,
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    Info<< "Reading initAlphaFieldDict" << endl;

    Foam::autoPtr<Foam::implicitFunction> func = implicitFunction::New
    (
           initAlphaFieldDict.get<word>("type"),
           initAlphaFieldDict
    );

    scalarField f(mesh.nPoints(),0.0);

    forAll(f,pi)
    {
        f[pi] = func->value(mesh.points()[pi]);
    };

    cutCellIso cutCell(mesh,f);

    DynamicList< List<point> > facePts;

    DynamicList <triSurface> surface;

    surfaceScalarField cellToCellDist(mag(mesh.delta()));

    forAll(alpha1,cellI)
    {
        label cellStatus = cutCell.calcSubCell(cellI,0.0);

        if(cellStatus == -1)
        {
            alpha1[cellI] = 1;
        }
        else if(cellStatus == 1)
        {
            alpha1[cellI] = 0;
        }
        else if(cellStatus == 0)
        {
            if(mag(cutCell.faceArea()) != 0)
            {
                alpha1[cellI]= max(min(cutCell.VolumeOfFluid(),1),0);
                if (writeVTK  && (mag(cutCell.faceArea()) >= 1e-14))
                {
                    facePts.append(cutCell.facePoints());
                }
            }
        }

    }

    if (writeVTK)
    {
        std::ostringstream os ;
        os << "AlphaInit" ;
        isoFacesToFile(facePts, os.str() , "AlphaInit");
    }

    ISstream::defaultPrecision(18);

    if(invert)
    {
        alpha1 = scalar(1) - alpha1;
    }

    alpha1.write();

    const scalarField& alpha = alpha1.internalField();
    Info << "sum(alpha*V) = " << gSum(mesh.V()*alpha)
     << ", 1-max(alpha1) = " << 1 - gMax(alpha)
     << "\t min(alpha1) = " << gMin(alpha) << endl;

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
