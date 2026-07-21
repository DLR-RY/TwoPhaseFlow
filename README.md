# twophaseflow

the TwoPhaseFlow Library adds new surface tension and phase change models to OpenFOAM and provides benchmark cases for verification

## Documentation

The available models and solvers are documentated in:

https://arxiv.org/abs/2103.00870

## Getting Started


### Prerequisites

Requires OpenFOAM v1812:

```
https://www.openfoam.com/download/release-history.php
```
The branch of2106 works with of2106

### Installing

```bash
    git clone https://github.com/DLR-RY/TwoPhaseFlow
    cd TwoPhaseFlow
    ./Allwmake
    ./get-gmsh.sh # will install gmsh version 306 as gmshv306
    # for AMR
    git submodule update --init --recursive
    cd modules/multiDimAMR/
    ./Allwmake
```
### running testsuite

make sure that the desired openfoam installation is sourced e.g. v1812 and that 
python is installed with a version >= 3.6 (miniconda is a great option, but anaconda works as well)

```bash
    python -m venv env # creats virtual python enviroments (optional step)
    pip install oftest

    py.test # runs the tests
    py.test --writeNSteps=1 run/ # test all testcases in run
```

## Authors

* **Henning Scheufler**

### adaptive mesh refinement with multiple regions

AMR with multiple regions does not work in version of1812 but it is fixed in newer versions.


To fix this apply the patch (assumes openfoam is already source):

```bash
    cp  patches/multiRegionAMR.patch $WM_PROJECT_DIR
    cp  patches/tableBase.patch $WM_PROJECT_DIR
    cp  patches/surfaceFieldValue.patch $WM_PROJECT_DIR
    cd $WM_PROJECT_DIR
    git apply multiRegionAMR.patch
    git apply tableBase.patch
    git apply surfaceFieldValue.patch

```
details see:

https://develop.openfoam.com/Development/openfoam/-/issues/1676

https://develop.openfoam.com/Development/openfoam/-/issues/1753
## License

This project is licensed under the GPL v3 License - see the [LICENSE.md](LICENSE.md) file for details



## running Benchmarks

```bash
    ./get-gmsh.sh # install gmsh
    pip install casefoam

```

The run/benchmark cases are run with


```bash
    cd run/benchmark/phaseChange/suckingInterface/
    python genCases.py # generates the study based and template case (here StefanProblem)
    ./Allrun # runs all the created testcases
    python getResults.py # to see results
```

Alternatively, the runAll.sh can be executed in the folder.

Note:

that some cases use the slurm queuing system and call `sbatch Allrun_Slurm` in the Allrun script, so you might need to modify it in the template case

