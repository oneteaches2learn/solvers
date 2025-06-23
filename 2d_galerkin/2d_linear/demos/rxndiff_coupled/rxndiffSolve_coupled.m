% RXNDIFFSOLVE_COUPLED
clear all; x = sym('x',[1 2],'real'); syms t; syms u; syms v;
% USER INPUTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% DOMAIN DATA
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


% TIME STEPPING DATA
T  = 0.1;
dt = 0.001;
eq.atEq = "break";
eq.tolerance = 1e-5;
dt = 0.001;


% PDE DATA
% specify coefficients
c1 = sin(x(1))^2 + t^2 + 1;
k  = cos(x(2))^2 + t^2 + 1;
r  = -10 * (u - v)^2;

% specify source
f = 0;

% specify BCs
bTypes = 'DDDD';
bTypes2 = 'D';

% specify Dirichlet conditions
u_D = 0;

% specify Neumann conditions
u_N = 0;

% specify Robin conditions
alpha = 1;
u_R = 1;

% specify initial condition
u_o = 0;


% ODE INFORMATION
g = 0;
v_o = 0;
c2 = 1;
s = @(u,v)(1 + u); % <~~~ NOTE: This is a placeholder value
s_const = 1;
constraints.vLower = NaN;
constraints.vUpper = NaN;
order = 1;


% RUN TRIAL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SETUP TRIAL
fprintf('Heat Trial Begun\n')

% assemble domain
fprintf(' Contructing Domain:'), tic
	%inc = Inclusion2d_circle(xLim_Y,yLim_Y,incRatio);
	%inc = Inclusion2d_square(xLim_Y,yLim_Y,incRatio);
	%dom = GalerkinAssembler2d_heat.assembleDomainGeometry(xLim_dom,yLim_dom,inc,eps);
	dom = GalerkinAssembler2d_rxndiff.assembleDomainGeometry(xLim_dom,yLim_dom);
	dom = GalerkinAssembler2d_rxndiff.assembleTimeStepping(dom,T,dt,eq);
	dom = dom.add_yline;
	%dom = dom.inclusionsON;
	dom = GalerkinAssembler2d_rxndiff.assembleBoundary(dom,bTypes,u_D,u_N,alpha,u_R,bTypes2); 
	dom = GalerkinAssembler2d_rxndiff.assembleMesh(dom,p,base); 
executionTime = toc; 
fprintf(' %f s\n',executionTime)

% assemble other parameters
fprintf(' Assembling Parameters:'), tic
	auxfun = GalerkinAssembler2d_rxndiff.assembleCoefficients(c1,k,r,f,u_o);
executionTime = toc; 
fprintf(' %f s\n',executionTime)

% Assemble ODE Data
fprintf(' Assembling ODE Data:'), tic
	data.g = g;
	data.vInit = v_o;
	data.cofs.c = c2;
	data.cofs.s = s;
	data.time.T = T;
	data.time.dt = dt;
	data.constraints = constraints;
	options.order = 1;
	ode = ODE(data,options);
executionTime = toc;
fprintf(' %f s\n',executionTime)

% RUN TRIAL
fprintf(' Generating Solution:'), tic
	prob = CoupledNewtonSolver2d_rxndiff(dom,auxfun,ode);
executionTime = toc; 
fprintf(' %f s\n',executionTime)