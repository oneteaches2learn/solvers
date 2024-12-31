% rxndiffSolve_demo 
clear all; x = sym('x',[1 2],'real'); syms t; syms u;
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
p  = 5;

% time stepping
T  = .1;
dt = 0.001;

% PDE INFORMATION
% specify coefficients
c = 1;
k = 1;
r = 100*u;
f = 100;

% specify boundary conditions
bcTypes_interior = 'DDDD';
bcTypes_exterior = 'R';
BCs = {0,0,0,0,{10,1}};

% specify initial condition
u_o = 0;

% ODE INFORMATION
f = 1;
v_o = 0;
s = 10;
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
	bcTypes = [bcTypes_interior, bcTypes_exterior];
	dom = GalerkinAssembler2d_rxndiff.assembleBoundary(dom_geo,bcTypes,BCs); 
	dom = GalerkinAssembler2d_rxndiff.assembleMesh(dom,p,base);
	dom = GalerkinAssembler2d_rxndiff.assembleTimeStepping(dom,T,dt);
	auxfun = GalerkinAssembler2d_rxndiff.assembleCoefficients(c,k,r,f,u_o);
executionTime = toc;
fprintf(' %f s\n',executionTime)

% Assemble ODE Data
fprintf(' Assembling ODE Data:'), tic
	data.f = f;
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

