%  Matlab mesh
% tutorial2, Created by Gmsh
% ASCII
clear msh;
msh.nbNod = 4;
msh.POS = [
0 0 0;
0.1 0 0;
0.1 0.3 0;
0 0.3 0;
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
 2 1 4 0
 2 4 3 0
];
msh.PNT =[
 1 0
 2 0
 3 0
 4 0
];
