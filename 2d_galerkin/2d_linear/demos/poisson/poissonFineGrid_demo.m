% PoissonMMS_demo
clear all; x = sym('x',[1 2],'real'); syms t;
% USER INPUTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% domain root
mshRoot = '/home/tyler/Software/MATLAB/solvers/2d_galerkin/2d_linear/gmsh_meshes/meshes/regular_mesh';

% fine grid parameters
pmin = 2;
pmax = 5;

% specify coefficients
k = 2 + sin(x(1) + x(2));
r = 2 + sin(x(1) + x(2));

% specify source
f = 10;

% specify BCs
bTypes_outer = 'RNDN';
bTypes_inner = 'R';

% specify Dirichlet conditions
u_D = 2 + sin(x(1) + x(2));

% specify Neumann conditions
u_N = 1;

% specify Robin conditions
beta = 1;
u_R = 1;



% FINE GRID TEST %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Poisson Fine Grid Test Begun\n')

% fine grid parameters
base = 2;
mmsparams = MMSParams(base,demo=0,timeOffset=4,timeFactor=2,pmin=pmin,pmax=pmax, ...
				meshInclusions='off',effectiveRegion='Omega)eps');
mmsparams.mshRoot = mshRoot;

% assemble domain
fprintf(' Contructing Domain:'), tic
	dom = Domain2d([0 1],[0 1]);
	dom = GalerkinAssembler2d_poisson.assembleBoundary(dom,bTypes_outer,u_D,u_N,beta,u_R,bTypes_inner); 
executionTime = toc; 
fprintf(' %f s\n',executionTime)

% assemble other parameters
fprintf(' Assembling Parameters:'), tic
	auxfun = GalerkinAssembler2d_poisson.assembleCoefficients(k,r,f);
executionTime = toc; 
fprintf(' %f s\n',executionTime)

% run trial
fprintf(' Generating Solution:'), tic
	fine = GalerkinFineGrid2d_poisson(dom,auxfun,mmsparams);
executionTime = toc; 
fprintf(' %f s\n',executionTime)

fine