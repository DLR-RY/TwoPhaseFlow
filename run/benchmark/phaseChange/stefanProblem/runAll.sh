#!/bin/sh
cd ${0%/*} || exit 1    # Run from this directory

python genCases.python


./Allrun