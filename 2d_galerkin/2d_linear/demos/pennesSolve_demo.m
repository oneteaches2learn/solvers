% pennesSolve_demo
clear all; x = sym('x',[1 2],'real'); syms t;
% USER INPUTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% domain bounds
xLim_dom = [0 1];
yLim_dom = [0 1];

% Y bounds
xLim_Y = [0 1];
yLim_Y = [0 1];

% number of inclusions
eps = 1/2;
incRatio = pi/2; % <~~~ incRatio = |delta Q| / |Y|

% mesh parameters
p = 4;
base = 2;

% specify coefficients
c = 1;
k = 1;
r = 1;
uStar = 1;

% specify source
f = 1;

% specify BCs
bTypes = 'DDDD';
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
	%inc = Inclusion2d_circle(xLim_Y,yLim_Y,incRatio);
	inc = Inclusion2d_square(xLim_Y,yLim_Y,incRatio);
	dom = GalerkinAssembler2d_pennes.assembleDomainGeometry(xLim_dom,yLim_dom,inc,eps);
	dom = GalerkinAssembler2d_pennes.assembleTimeStepping(dom,T,dt,eq);
	dom = dom.add_yline;
	dom = dom.inclusionsOFF;
	dom = GalerkinAssembler2d_pennes.assembleBoundary(dom,bTypes,u_D,u_N,beta,u_R,bTypes2); 
	dom = GalerkinAssembler2d_pennes.assembleMesh(dom,p,base); 
executionTime = toc; 
fprintf(' %f s\n',executionTime)

% assemble other parameters
fprintf(' Assembling Parameters:'), tic
	auxfun = GalerkinAssembler2d_pennes.assembleCoefficients(c,k,r,uStar,f,u_o);
executionTime = toc; 
fprintf(' %f s\n',executionTime)

% RUN TRIAL
fprintf(' Generating Solution:'), tic
	prob = GalerkinSolver2d_pennes(dom,auxfun);
executionTime = toc; 
fprintf(' %f s\n',executionTime)