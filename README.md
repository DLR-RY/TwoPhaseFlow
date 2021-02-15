# twophaseflow

the twophaseflow adds new surface tension and phase change models to OpenFOAM and provides benchmark cases for verification

## Getting Started


### Prerequisites

Requires OpenFOAM v1812:

```
https://www.openfoam.com/download/release-history.php
```

### Installing

```
git clone https://github.com/DLR-RY/TwoPhaseFlow
cd twophaseflow
./Allwmake
# for AMR
git submodule update --init
cd modules/multiDimAMR/
./Allwmake
```

## Authors

* **Henning Scheufler**

### adaptive mesh refinement with multiple regions

AMR with multiple regions does not work in version of1812 but it is fixed in newer versions.


To fix this apply the patch (assumes openfoam is already source):

```
cp  patches/multiRegionAMR.patch $WM_PROJECT_DIR
cd $WM_PROJECT_DIR
git apply multiRegionAMR.patch

```
details see:

https://develop.openfoam.com/Development/openfoam/-/issues/1676
## License

This project is licensed under the GPL v3 License - see the [LICENSE.md](LICENSE.md) file for details




