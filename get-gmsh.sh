#!/bin/sh
cd "${0%/*}" || exit                                # Run from this directory

wget https://gmsh.info/bin/Linux/gmsh-3.0.6-Linux64.tgz
tar zxvf gmsh-3.0.6-Linux64.tgz
cp gmsh-3.0.6-Linux64/bin/gmsh $FOAM_USER_APPBIN/gmshv306
rm gmsh-3.0.6-Linux64.tgz
