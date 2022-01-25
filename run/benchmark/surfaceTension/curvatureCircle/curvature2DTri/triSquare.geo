nx=10;

Point(1) = {-0, -0.05, 0, 1e+22};
  Point(2) = {2, -0.05, 0, 1e+22};
  Point(3) = {2, -0.05, 2, 1e+22};
  Point(4) = {0, -0.05, 2, 1e+22};
  Line(1) = {1, 2};
  Line(2) = {2, 3};
  Line(3) = {3, 4};
  Line(4) = {4, 1};
  Line Loop(6) = {4, 1, 2, 3};

  Plane Surface(6) = {6};
  Physical Volume("internal") = {1};

 Transfinite Line{1} = nx;
 Transfinite Line{2} = nx;
 Transfinite Line{3} = nx;
 Transfinite Line{4} = nx;

  Extrude {0, 0.1, 0} {
   Surface{6};
   Layers{1};
   Recombine;
  }
//  Physical Surface("front") = {28};
//  Physical Surface("back") = {6};
//  Physical Surface("bottom") = {27};
//  Physical Surface("left") = {15};
//  Physical Surface("top") = {19};
//  Physical Surface("right") = {23};

  Physical Surface("back") = {6};
  Physical Surface("rim") = {27, 15, 19, 23};
  Physical Surface("front") = {28};
