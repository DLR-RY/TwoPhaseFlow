#!/bin/sh
cd ${0%/*} || exit 1    # Run from this directory

python genCases.py

cat Allrun | tr -d "&" > tmp
mv tmp Allrun
chmod u+x Allrun

./Allrun

python getResults.py