#!/bin/sh
cd "${0%/*}" || exit                                # Run from this directory

wget https://gmsh.info/bin/Linux/gmsh-4.9.3-Linux64.tgz
tar zxvf gmsh-4.9.3-Linux64.tgz
cp gmsh-4.9.3-Linux64/bin/gmsh $FOAM_USER_APPBIN/gmshv493
rm gmsh-4.9.3-Linux64.tgz
