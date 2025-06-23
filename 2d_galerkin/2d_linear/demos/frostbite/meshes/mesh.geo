 // Inputs
squareSide = 1; 
meshThickness = squareSide; 
gridsize = squareSide * 1000;

// Geometry
Point(1) = {0, 0, 0, gridsize};
Point(2) = {squareSide, 0, 0, gridsize};
Point(3) = {squareSide, squareSide, 0, gridsize};
Point(4) = {0, squareSide, 0, gridsize};
Line(1) = {1, 2};				// bottom line
Line(2) = {2, 3};				// right line
Line(3) = {3, 4};				// top line
Line(4) = {4, 1};				// left line
Line Loop(5) = {1, 2, 3, 4}; 	
Plane Surface(6) = {5};