%  Matlab mesh
% mesh, Created by Gmsh
% ASCII
clear msh;
msh.nbNod = 5;
msh.POS = [
0 0 0;
1 0 0;
1 1 0;
0 1 0;
0.5 0.5 0;
];
msh.MAX = max(msh.POS);
msh.MIN = min(msh.POS);
msh.LINES =[
 1 2 0
 2 3 0
 3 4 0
 4 1 0
];
msh.TRIANGLES =[
 1 2 5 0
 4 1 5 0
 2 3 5 0
 3 4 5 0
];
msh.PNT =[
 1 0
 2 0
 3 0
 4 0
];
