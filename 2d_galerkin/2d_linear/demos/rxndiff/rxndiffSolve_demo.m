% HEATSOLVE_DEMO
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
p = 5;
base = 2;

% specify coefficients
c = 1;
k = 1;

% specify source
f = 2;

% specify BCs
bTypes = 'RRRR';
bTypes2 = 'R';

% specify Dirichlet conditions
u_D = 0;

% specify Neumann conditions
u_N = 0;

% specify Robin conditions
alpha = cos(pi/2 * x(1)) * cos(pi/2 * x(2));
u_R = 1;

% specify Dynamic conditions
beta = 1;
gamma = 1;
u_T = 1;

% time-stepping parameters
u_o = 0;
T  = 1;
dt = 0.001;
eq.atEq = "break";
eq.tolerance = 1e-5;
dt = 0.001;


% RUN TRIAL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SETUP TRIAL
fprintf('Heat Trial Begun\n')

% assemble domain
fprintf(' Contructing Domain:'), tic
	inc = Inclusion2d_circle(xLim_Y,yLim_Y,incRatio);
	%inc = Inclusion2d_square(xLim_Y,yLim_Y,incRatio);
	dom = GalerkinAssembler2d_heat.assembleDomainGeometry(xLim_dom,yLim_dom,inc,eps);
	dom = GalerkinAssembler2d_heat.assembleTimeStepping(dom,T,dt,eq);
	dom = dom.add_yline;
	%dom = dom.inclusionsON;
	dom = GalerkinAssembler2d_heat.assembleBoundary(dom,bTypes,u_D,u_N,alpha,u_R,bTypes2); 
	dom = GalerkinAssembler2d_heat.assembleMesh(dom,p,base); 
executionTime = toc; 
fprintf(' %f s\n',executionTime)

% assemble other parameters
fprintf(' Assembling Parameters:'), tic
	auxfun = GalerkinAssembler2d_heat.assembleCoefficients(c,k,0,f,u_o);
executionTime = toc; 
fprintf(' %f s\n',executionTime)

% RUN TRIAL
fprintf(' Generating Solution:'), tic
	prob = GalerkinSolver2d_heat(dom,auxfun);
executionTime = toc; 
fprintf(' %f s\n',executionTime)