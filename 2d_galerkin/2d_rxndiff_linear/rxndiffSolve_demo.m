% rxndiffSolve_demo 
clear all; x = sym('x',[1 2],'real'); syms t;
% USER INPUTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% domain bounds
xLim_dom = [0 3];
yLim_dom = [0 1];

% Y bounds
xLim_Y = [0 1/2];
yLim_Y = [0 1];

% number of inclusions
eps = 1/3;
incRatio = pi/2; % <~~~ incRatio = |delta Q| / |Y|

% mesh parameters
p = 5;
base = 2;

% specify coefficients
c = 1;
k = 1;
r = 1;

% specify source
f = 1;

% specify BCs
bTypes = {'D' 'D' 'D' 'D'};
bTypes2 = 'D';

% specify Dirichlet conditions
u_D = 0;

% specify Neumann conditions
u_N = 1;

% specify Robin conditions
beta = 1;
u_R = 1;

% time-stepping parameters
u_o = 0;
T  = 1;
dt = 0.001;
eq.atEq = "break";
eq.tolerance = 1e-5;


% RUN TRIAL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SETUP TRIAL
fprintf('RxnDiff Trial Begun\n')

% assemble domain
fprintf(' Contructing Domain:'), tic
	%inc = Inclusion2d_circle(xLim_Y,yLim_Y,incRatio);
	inc = Inclusion2d_square(xLim_Y,yLim_Y,incRatio);
	dom = GalerkinRxndiff2d_assembler.assembleDomainGeometry(xLim_dom,yLim_dom,inc,eps);
	dom = GalerkinRxndiff2d_assembler.assembleBoundary(dom,bTypes,u_D,u_N,beta,u_R,bTypes2); 
	dom = GalerkinRxndiff2d_assembler.assembleMesh(dom,p,base); 
executionTime = toc; 
fprintf(' %f s\n',executionTime)

% assemble other parameters
fprintf(' Assembling Parameters:'), tic
	cofs   = GalerkinRxndiff2d_assembler.assembleCoefficients(c,k,r);
	uInit  = GalerkinRxndiff2d_assembler.assembleInitialCondition(u_o);
	time   = GalerkinRxndiff2d_assembler.assembleTimeStepping(T,dt);
	source = GalerkinRxndiff2d_assembler.assembleSource(f);
executionTime = toc; 
fprintf(' %f s\n',executionTime)

% RUN TRIAL
fprintf(' Generating Solution:'), tic
	prob = GalerkinRxndiff2d_solver(dom,time,cofs,uInit,source);
executionTime = toc; 
fprintf(' %f s\n',executionTime)
