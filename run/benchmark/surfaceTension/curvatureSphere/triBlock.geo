nx=replaceNx;
nz=replaceNz;


// Gmsh project created on Mon Jun 12 15:29:01 2017
SetFactory("OpenCASCADE");
//+
Block(1) = {0, 0, 0, 1, 1, 1};
//+
Transfinite Line {11, 9, 12, 10} = nx Using Progression 1;
//+
Transfinite Line {7, 6, 5, 8, 3, 4, 1, 2} = nx Using Progression 1;
//+
Physical Surface("wall") = {6, 2, 4, 3, 5, 1};
Physical Volume("internal") = {1};

