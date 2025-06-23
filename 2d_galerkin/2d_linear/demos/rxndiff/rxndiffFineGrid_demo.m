% rxndiffMMS_demo
clear all; x = sym('x',[1 2],'real'); syms t;
% USER INPUTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% domain root
mshRoot = '/home/tyler/Software/MATLAB/solvers/2d_galerkin/2d_linear/gmsh_meshes/meshes/regular_mesh';

% fine grid parameters
pmin = 2;
pmax = 5;


% specify coefficients
%c = 2 + sin(x(1) + x(2)) * t + t^2; 
%k = 2 + sin(x(1) + x(2)) * t + t^2; 
c = 1;
k = 1;
r = 1;

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

% specify time parameters
T = 1;
uInit = u_D;


% FINE GRID TEST %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Reaction Diffusion Fine Grid Test Begun\n')

% assemble inputs
base = 2;
mmsparams = MMSParams(base,demo=0,timeOffset=pmin,timeFactor=2,pmin=pmin,pmax=pmax, ...
                meshInclusions='off',effectiveRegion='Omega)eps');
mmsparams.mshRoot = mshRoot;
mmsparams.pfine_add = 3;

% assemble domain
fprintf(' Contructing Domain:'), tic
    dom = Domain2d([0 1],[0 1]);
    dom = GalerkinAssembler2d_rxndiff.assembleBoundary(dom,bTypes_outer,u_D,u_N,beta,u_R,bTypes_inner); 
    dom.time = TimeStepping(T,1);
executionTime = toc; 
fprintf(' %f s\n',executionTime)

% assemble other parameters
fprintf(' Assembling Parameters:'), tic
    auxfun = GalerkinAssembler2d_rxndiff.assembleCoefficients(c,k,r,f,uInit);
executionTime = toc; 
fprintf(' %f s\n',executionTime)

% run trial
fprintf(' Generating Solution:'), tic
    fine = GalerkinFineGrid2d_rxndiff(dom,auxfun,mmsparams);
executionTime = toc; 
fprintf(' %f s\n',executionTime)

fine