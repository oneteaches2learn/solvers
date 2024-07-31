% DEMO3_DRIVER 
% Generates a manufactured solution and the files necessary to replicate that
% manufactured solution using the 51 Lines of Matlab code. The purpose is to
% compare the vectors and tensors that my code is producing against the ACF
% code. 

clear all; x = sym('x',[1 2],'real'); syms t;
% USER INPUTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% domain bounds
xLim = [0 1];
yLim = [0 1];

% time stepping
T = 1;

% number of inclusions
N_x = 3;
N_y = 3;
alpha = 1; % <~~~ alpha = |delta Q| / |Y|

% generate mesh
q = 5;
base = 2;
h = base^-q;

% specify BCs
bTypes = {'D' 'D' 'D' 'D'};
bTypes2 = 'D';

% specify coefficients
p = 1;
k = 1;
r = 0;

% specify desired result
uTrue = sin(pi*x(1))*sin(pi*x(2))*t;

mms_test = 0;


% MMS TEST %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('MMS Test Begun\n')

% assemble inputs
bound  = PuncturedBoundary2d(bTypes,{@()(0.0),@()(0.0),@()(0.0),@()(0.0)},bTypes2,{@()(0.0)});
auxfun = HeatAuxFunctions2d(p,k,r,@()(0.0),uTrue);
time   = TimeStepping(T,1);

% build domain
fprintf('Initialization\n')
fprintf(' Contructing Domain:'), tic
incMod = InclusionModule1(alpha);
%dom    = PuncturedDomain2d(xLim,yLim,N_x,N_y,incMod);
dom    = Domain2d(xLim,yLim);
dom = dom.setEdgeBCTypes(bound);
executionTime = toc; 
fprintf(' %f s\n',executionTime)


if 1 == 1
% run mms test
if mms_test == 1
	mms = GalerkinHeat2d_mms(dom,time,auxfun,errType="Linfty(L2)")

% run demo test
elseif mms_test == 0
	mms = GalerkinHeat2d_mms(dom,time,auxfun,demo=q);
	prob = mms.problems{1};
end

end


% RUN ACF CODE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
path = "/home/tyler/Software/MATLAB/solvers/fem2d_heat/";

% generate node coordinates
coordinates = zeros(prob.domain.nNodes,3);
coordinates(:,1) = [1:1:prob.domain.nNodes];
coordinates(:,2:3) = prob.domain.Mesh.Nodes';
save(path + "coordinates.dat",'coordinates','-ascii','-tabs');

% generate dirichlet segments
dirichlet = [];
curRow = 1;
for i = 1:prob.domain.NumEdges
	if prob.domain.edges(i).boundaryType == 'D'
		for j = 1:prob.domain.edges(i).nNodes-1
			dirichlet(curRow,1) = curRow;
			dirichlet(curRow,2) = prob.domain.edges(i).nodes(j);
			dirichlet(curRow,3) = prob.domain.edges(i).nodes(j+1);
			curRow = curRow + 1;
		end
	end
end
save(path + "dirichlet.dat",'dirichlet','-ascii','-tabs');

% generate neumann segments
neumann = [];
curRow = 1;
for i = 1:prob.domain.NumEdges
	if prob.domain.edges(i).boundaryType == 'N'
		for j = 1:prob.domain.edges(i).nNodes-1
			neumann(curRow,1) = curRow;
			neumann(curRow,2) = prob.domain.edges(i).nodes(j);
			neumann(curRow,3) = prob.domain.edges(i).nodes(j+1);
			curRow = curRow + 1;
		end
	end
end
save(path + "neuamnn.dat",'neumann','-ascii','-tabs');

% generate elements3
elements3 = zeros(size(prob.domain.Mesh.Elements,2),3);
elements3(:,1) = [1:1:size(prob.domain.Mesh.Elements,2)];
elements3(:,2:4) = prob.domain.Mesh.Elements';
save(path + "elements3.dat",'elements3','-ascii','-tabs');

% generate f
func = prob.f;
fileID = fopen(path + "f.m",'w');
fprintf(fileID,"function VolumeForce = f(x,t)\n\n");
fprintf(fileID,"\tx1 = x(1);\n\tx2 = x(2);\n");
fprintf(fileID,"\tfunc = " + func2str(func) + ";\n");
fprintf(fileID,"\tVolumeForce = func(x1,x2,t);\n\n");
fprintf(fileID,"end");


% generate u_d
func = prob.domain.edges(1).boundaryCondition;
fileID = fopen(path + "u_d.m",'w');
fprintf(fileID,"function DirichletBoundaryValue = u_d(x,t)\n\n");
fprintf(fileID,"\tx1 = x(1);\n\tx2 = x(2);\n\t");
fprintf(fileID,"\tfunc = " + func2str(func) + ";\n\n");
fprintf(fileID,"\tDirichletBoundaryValue = func(x1,x2,t);\n\n");
fprintf(fileID,"end");

fem2d_heat

