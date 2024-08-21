% PENNESSOLVE_DEMO

clear all; x = sym('x',[1 2],'real'); syms t;
% USER INPUTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% domain bounds
xLim = [0 1];
yLim = [0 1];

% number of inclusions
N_x = 2;
N_y = N_x;

% mesh parameters
p = 4;
base = 2;
incRatio = 2; % <~~~ incRatio = |delta Q| / |Y|
eps = 1;

% specify coefficients
c = 1;
k = 1;
r = 1;
uStar = 1;

% specify source
f = 1;

% specify BCs
bTypes = {'D' 'D' 'D' 'D'};
bTypes2 = 'R';

% specify Dirichlet conditions
u_D = 0;

% specify Neumann conditions
u_N = 1;

% specify Robin conditions
beta = 1;
u_R = uStar;

% time-stepping parameters
u_o = 0;
T  = 1;
dt = 0.001;
eq.atEq = "break";
eq.tolerance = 1e-5;


% RUN TRIAL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SETUP TRIAL
fprintf('Pennes Trial Begun\n')

% assemble domain
fprintf(' Contructing Domain:'), tic
	incMod = InclusionModule_square(incRatio,eps);
	dom = GalerkinPennes2d_assembler.assembleDomainGeometry(xLim,yLim,incMod);
	dom = GalerkinPennes2d_assembler.assembleBoundary(dom,bTypes,u_D,u_N,beta,u_R,bTypes2); 
	dom = GalerkinPennes2d_assembler.assembleMesh(dom,p,base); 
executionTime = toc; 
fprintf(' %f s\n',executionTime)

% assemble other parameters
fprintf(' Assembling Parameters:'), tic
	cofs   = GalerkinPennes2d_assembler.assembleCoefficients(c,k,r,uStar);
	uInit  = GalerkinPennes2d_assembler.assembleInitialCondition(u_o);
	time   = GalerkinPennes2d_assembler.assembleTimeStepping(T,dt,eq);
	source = GalerkinPennes2d_assembler.assembleSource(f);
executionTime = toc; 
fprintf(' %f s\n',executionTime)

% RUN TRIAL
fprintf(' Generating Solution:'), tic
	prob = GalerkinPennes2d_solver(dom,time,cofs,uInit,source);
executionTime = toc; 
fprintf(' %f s\n',executionTime)
