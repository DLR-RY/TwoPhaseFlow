#!/bin/sh
cd "${0%/*}" || exit                                # Run from this directory

git submodule update --init --recursive

cd modules/multiDimAMR
    ./Allwmake
cd ../..

