import sys
import subprocess

# Picture to be meshed
filename = str(sys.argv[1])
filename = f'"{filename}"'

# Write GEO file as input for Gmsh
header = f'If(PostProcessing.NbViews == 0) \n \
\t Merge {filename}; \n\
\t Plugin(ModifyComponents).Expression0 = "1 + v0^3 * 10"; \n\
\t Plugin(ModifyComponents).Run; \n\
EndIf\n'

body = """
Background Mesh View[0];
w = View[0].MaxX;
h = View[0].MaxY;

Point(1)={0,0,0,w};
Point(2)={w,0,0,w};
Point(3)={w,h,0,w};
Point(4)={0,h,0,w};
Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};
Line Loop(5) = {3,4,1,2};
Plane Surface(6) = {5};

DefineConstant[
  algo = {Mesh.Algorithm, AutoCheck 0, GmshOption "Mesh.Algorithm",
    Choices{1="MeshAdapt", 2="Automatic", 5="Delaunay", 6="Frontal-Delaunay", 
    7="BAMG", 8="Frontal-Delaunay for Quads", 9="Packing of Parallelograms",
    11="Quasi-structured Quad"},
    Name "Meshing parameters/Algorithm"},

  sizeFact = {Mesh.CharacteristicLengthFactor, AutoCheck 0,
    GmshOption "Mesh.CharacteristicLengthFactor", Min 0.1, Max 10, Step 0.1,
    Name "Meshing parameters/Element size factor"},

  sizeMin = {Mesh.CharacteristicLengthMin, AutoCheck 0,
    GmshOption "Mesh.CharacteristicLengthMin", Min w/100, Max w, Step 0.1,
    Name "Meshing parameters/Minimum element size"},

  save = {StrCat("View.ShowScale=0;Print '", CurrentDirectory, "outputFile.png';"),
    AutoCheck 0, Macro "GmshParseString",
    Name "Save PNG"}
];

Mesh.ColorCarousel = 0;
Mesh.Color.Triangles = Black;
Mesh.Color.Quadrangles = Black;
Mesh.RecombineAll = (algo == 8);
Solver.AutoMesh = 2;"""

geoFile = header + body

def write_input():
    f = open('input.geo', "w")
    f.write(geoFile)
    f.close()

def run_mesh_gui():
    subprocess.run(["gmsh", "input.geo"])
    subprocess.run(["gmsh", "input.msh"])

# Run
write_input()
run_mesh_gui()

