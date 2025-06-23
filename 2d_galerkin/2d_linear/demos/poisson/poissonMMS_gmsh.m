% PoissonMMS_demo
clear all; x = sym('x',[1 2],'real'); syms t;
% USER INPUTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% domain root
mshRoot = '/home/tyler/Software/MATLAB/solvers/2d_galerkin/2d_linear/gmsh_meshes/meshes/irregular_mesh';

% time stepping
T = 1;

% mms parameters
demo = 0;

% specify BCs
bTypes_outer = 'DDDD';
bTypes_inner = 'R';

% specify coefficients
k = 2 + sin(x(1) + x(2));
r = 2 + sin(x(1) + x(2));

% specify desired result
%uTrue = sin(pi/2 * x(1)) * sin(pi/2 * x(2));
%uTrue = sin(pi * x(1)) * sin(pi * x(2));
%uTrue = cos(2 * pi * x(1)) * cos(2 * pi * x(2));
uTrue = 1;
uTrue = sin(2 * pi * x(1)) * sin(2 * pi * x(2));

% MMS TEST %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('MMS Test Begun\n')

% assemble inputs
auxfun    = ManufacturedFunctions2d_poisson(k,r,uTrue);
base = 2;
mmsparams = MMSParams(base,demo=demo,timeOffset=4,timeFactor=2,pmin=4,pmax=8, ...
				meshInclusions='off',effectiveRegion='Omega_eps');
mmsparams.mshRoot = mshRoot;

% build dummy domain as container for boundary conditions
fprintf('Initialization\n')
fprintf(' Contructing Domain:'), tic
	xLim_dom = [0 1];
	yLim_dom = [0 1];
	dom = Domain2d(xLim_dom,yLim_dom);
	dom = dom.setBCTypes([bTypes_outer,bTypes_inner]);
executionTime = toc; 
fprintf(' %f s\n',executionTime)


% run mms test
if demo == 0
	mms = GalerkinMMS2d_poisson(dom,auxfun,mmsparams)

% run demo test
else
	mms = GalerkinMMS2d_poisson(dom,auxfun,mmsparams);
	prob = mms.problems{1};
end