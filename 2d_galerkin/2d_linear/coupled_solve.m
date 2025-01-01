% rxndiffSolve_demo 
clear all; x = sym('x',[1 2],'real'); syms t; syms u; syms v;
% USER INPUTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% DOMAIN INFORMATION
% domain bounds
xLim_dom = [0 1];
yLim_dom = [0 1];

% Y bounds
xLim_Y = [0 1];
yLim_Y = [0 1];

% number of inclusions
eps = 1/4;
incRatio = pi/2; % <~~~ incRatio = |delta Q| / |Y|

% domain parameters
base = 2;
p  = 3;

% time stepping
T  = 3;
dt = 0.05;

% PDE INFORMATION
% specify coefficients
c = 1;
k = 1;
r = u - v;
r = 0;
f = 1;

% manufacture solution for testing
%uTrue = sin(pi / 2 * x(1)) * sin(pi / 2 * x(2)) * (1 - exp(-10 * t));
%uTrue = sin(pi / 2 * x(1)) * sin(pi / 2 * x(2));
%uTrue = t;
%uTrue_t = diff(uTrue,t);
%uTrue_laplacian = divergence(gradient(uTrue,x),x);
%f = uTrue_t - uTrue_laplacian + (uTrue - v);

% specify boundary conditions
bcTypes_exterior = 'DDDD';
bcTypes_interior = 'R';
%bc1 = sin(pi / 2 * x(1)) * sin(pi / 2 * 0);
%bc2 = sin(pi / 2 * 1) * sin(pi / 2 * x(2));
%bc3 = sin(pi / 2 * x(1)) * sin(pi / 2 * 1);
%bc4 = sin(pi / 2 * 0) * sin(pi / 2 * x(2));
bc1 = 1;
bc2 = 1;
bc3 = 1;
bc4 = 1;
%bc1 = t;
%bc2 = t;
%bc3 = t;
%bc4 = t;
BCs = {bc1,bc2,bc3,bc4,{10,1}};

% specify initial condition
%u_o = sin(pi / 2 * x(1)) * sin(pi / 2 * x(2));
u_o = 1;


% ODE INFORMATION
g = 0;
v_o = 0;
s = 0;
order = 1;


% BLACK BOX %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SETUP TRIAL
fprintf('Solving Nonlinear RxnDiff\n')

% Construct Domain Geometry
fprintf(' Contructing Domain Geometry:'), tic
	%inc = Inclusion2d_square(xLim_Y,yLim_Y,incRatio);
	%inc = Inclusion2d_circle(xLim_Y,yLim_Y,incRatio);
	%dom_geo = GalerkinAssembler2d_rxndiff.assembleDomainGeometry(xLim_dom,yLim_dom,inc,eps);
	dom_geo = GalerkinAssembler2d_rxndiff.assembleDomainGeometry(xLim_dom,yLim_dom);
executionTime = toc; 
fprintf(' %f s\n',executionTime)

% Assemble PDE Data
fprintf(' Assembling PDE Data:'), tic
	bcTypes = [bcTypes_exterior, bcTypes_interior];
	dom = GalerkinAssembler2d_rxndiff.assembleBoundary(dom_geo,bcTypes,BCs); 
	dom = GalerkinAssembler2d_rxndiff.assembleMesh(dom,p,base);
	dom = GalerkinAssembler2d_rxndiff.assembleTimeStepping(dom,T,dt);
	auxfun = GalerkinAssembler2d_rxndiff.assembleCoefficients(c,k,r,f,u_o);
executionTime = toc;
fprintf(' %f s\n',executionTime)

% Assemble ODE Data
fprintf(' Assembling ODE Data:'), tic
	data.g = g;
	data.vInit = v_o;
	data.cofs.s = s;
	data.time.T = T;
	data.time.dt = dt;
	options.order = order;
	ode = ODE(data,options);
executionTime = toc;
fprintf(' %f s\n',executionTime)


% Solve
fprintf(' Solving:'), tic
	%prob = NewtonGalerkinSolver2d_rxndiff(dom,auxfun);
	prob = CoupledNewtonSolver2d_rxndiff(dom,auxfun,ode);
	%prob = NewtonGalerkinSolver2d_poisson(dom,auxfun);
	%prob = NewtonGalerkinSolver2d_heat(dom,auxfun);
executionTime = toc; 
fprintf(' %f s\n',executionTime)
