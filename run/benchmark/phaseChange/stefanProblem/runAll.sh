#!/bin/bash
cd ${0%/*} || exit 1                        # Run from this directory

python genCases.python

while read fname ; do
  if [[ !("$fname" =~ "baseCase") ]]
  then
    ($fname) &
  fi
done < <(find Cases -name "Allrun")

# wait for all process to finish
wait

python getData.py
python plotData.py
