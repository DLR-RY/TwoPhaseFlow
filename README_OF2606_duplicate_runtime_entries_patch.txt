OF2606 runtime-selection correction
===================================

The previous archive incorrectly wrapped several runtime-selection registrations
in `#if OPENFOAM < 2606` guards. That prevents OpenFOAM v2606 from seeing
important schemes such as `plicRDF` at run time, causing errors like:

    Unknown reconstructionSchemes type plicRDF

This corrected archive registers these classes unconditionally:

    reconstructionSchemes: plicRDF, gradAlpha, isoAlpha
    sampledSurface: sampledInterface as interface

The lnInclude copies included in this zip were updated as well so the archive is
self-consistent before running wmake/wmakeLnInclude.

Recommended rebuild:

    ./Allwclean
    ./Allwmake 2>&1 | tee log.Allwmake
