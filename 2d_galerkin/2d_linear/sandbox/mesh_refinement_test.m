% mesh_refinement_test
%  test of MATLAB's inbuilt mesh refinement tools. In particular, I want to see
%  if the mesh can be refined while keeping all the nodes from the original
%  mesh.

%{
% make a mesh by calling the basic poisson solve demo
poissonSolve_demo

% get inputs for refinement tool
dl_mat = prob.domain.boundary.dl.mat;
[p,e,t] = meshToPet(prob.domain.mesh.Mesh);

% refine mesh
[p1,e1,t1] = refinemesh(dl_mat,p,e,t);

% make dummy model and generate FEmesh
model = createpde();
geometryFromMesh(model,p1,t1(1:3,:));
FEmsh = model.Mesh;

% make Mesh2d object
msh = Mesh2d(FEmsh);

% plot meshes
figure("visible","on")
subplot(2,1,1)
prob.domain.mesh.plot;
subplot(2,1,2)
msh.plot;
%}

%prob2 = prob;
%prob2.domain = prob2.domain.refineMesh;

%{
clear all
dom = Domain2d([0 1],[0 1]);

dom = dom.setMesh(2,2);

dom2 = dom.refineMesh(1);
dom3 = dom.refineMesh(2);
dom4 = dom.refineMesh(3);
%}

clear all
dom = Domain2d([0 1],[0 1]);

fineParms = struct('pmin',1,'pmax',3,'base',2,'timeOffset',0,'timeFactor',1,...
	'effectiveRegion','all','meshInclusions','off','mshRoot','');

fine = GalerkinFineGrid2d(dom,[],fineParms);