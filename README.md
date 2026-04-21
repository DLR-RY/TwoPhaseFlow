# TwoPhaseFlow

The TwoPhaseFlow library adds new surface tension and phase change models to
OpenFOAM and provides benchmark cases for verification.

## Documentation

The available models and solvers are documented in the paper:

Scheufler, H., & Roenby, J. (2023). "TwoPhaseFlow: A Framework for Developing
Two Phase Flow Solvers in OpenFOAM". OpenFOAM® Journal, 3, 200–224.
https://doi.org/10.51560/ofj.v3.80

For examples of applications and extensions of the TwoPhaseFlow library, have a
look at the citing articles on [google scholar](https://scholar.google.com/scholar?oi=bibs&hl=en&cites=5885558257260540925,815283593763527218&as_sdt=5).

## Getting Started


### Prerequisites

Requires an installed OpenCFD version of OpenFOAM (v1812 or later) from
https://www.openfoam.com.

The master branch currently (april 2026) compiles with OpenFOAM-v2406 to v2512.
For older OpenFOAM versions (down to v1812) please checkout the appropriate
branch. Don't expect older branches to be regularly updated or bug fixed. Use a
newer OpenFOAM version if possible.

OpenFOAM.org versions are not supported.

### Installing

```bash
    git clone https://github.com/DLR-RY/TwoPhaseFlow
    cd TwoPhaseFlow
    # To compile e.g. with OpenFOAM-v2112 checkout the appropriate branch with:
    # git checkout of2112
    ./Allwmake
    ./get-gmsh.sh # will install gmsh version 306 as gmshv306
    # For AMR with DLB:
    git submodule update --init --recursive
    cd modules/multiDimAMR/
    ./Allwmake
```
### Running testsuite

Make sure that the desired OpenFOAM installation is sourced e.g. v2512 and that 
python is installed with a version >= 3.6 (miniconda is a great option, but
anaconda works as well)

```bash
    python -m venv env # creats virtual python enviroments (optional step)
    pip install oftest

    py.test # runs the tests
    py.test --writeNSteps=1 run/ # test all testcases in run
```

## Authors

* **Henning Scheufler**

## License

This project is licensed under the GPL v3 License - see the
[LICENSE.md](LICENSE.md) file for details.


## Running benchmarks

```bash
    ./get-gmsh.sh # Install gmsh
    pip install casefoam

```

The run/benchmark cases are run with


```bash
    cd run/benchmark/phaseChange/suckingInterface/
    python genCases.py # Generates the study cases based and template case (here StefanProblem)
    ./Allrun # Runs all the created testcases
    python getResults.py # See results
```

Alternatively, the runAll.sh can be executed in the folder.

Note:
Some cases use the slurm queuing system and call `sbatch Allrun_Slurm` in the
Allrun script, so you might need to modify it in the template case.

## Known issues

### Duplicity entry warnings in log file

For some of the solvers, depending on the libraries loaded with libs(...) in the controlDict, the log file may contain:
```
Duplicate entry isoAlpha in runtime table reconstructionSchemes
[stack trace]
=============
#1      platforms/linux64GccDPInt32Opt/lib/libgeometricVoF.so
...
```
This is because both the TwoPhaseFlow library ``libVoF.so`` and the OpenFOAM library ``libgeometricVoF.so`` add the reconstrucionSchemes to the runTime selection table. If for instance some loaded library links ``libincompressibleMultiphaseSystems.so``, this loads ``libgeometricVoF.so`` and we get the conflict. The solver will, however, use the ``libVoF.so`` version as it should so the warnings can be safely ignored.

### Adaptive mesh refinement with multiple regions with OpenFOAM-v1812

AMR with multiple regions does not work in version of1812 but it is fixed in
newer versions.

To fix this, apply the patches (assumes OpenFOAM is already source):

```bash
    cp  patches/multiRegionAMR.patch $WM_PROJECT_DIR
    cp  patches/tableBase.patch $WM_PROJECT_DIR
    cp  patches/surfaceFieldValue.patch $WM_PROJECT_DIR
    cd $WM_PROJECT_DIR
    git apply multiRegionAMR.patch
    git apply tableBase.patch
    git apply surfaceFieldValue.patch
    ./Allwmake
```
Important note: These patches are only possible if you have write access to your
OpenFOAM-v1812 installation. After applying the pathes your OpenFOAM-v1812
installation will no longer be identical to the official OpenFOAM-v1812 release.

For details see:

https://gitlab.com/openfoam/core/openfoam/-/work_items/1676

https://gitlab.com/openfoam/core/openfoam/-/work_items/1753
