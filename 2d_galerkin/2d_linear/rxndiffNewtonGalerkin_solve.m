% rxndiffSolve_demo 
%clear all; 
x = sym('x',[1 2],'real'); syms t; syms u;
% USER INPUTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
T  = .1;
dt = 0.001;

% specify coefficients
c = 10;
k = 1;
r = 0;
f = 10;

% specify boundary conditions
bcTypes_interior = 'DDDD';
bcTypes_exterior = 'R';
BCs = {0,0,0,0,{10,1}};

% specify initial condition
u_o = 0;


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

% Assemble Data
fprintf(' Assembling Data:'), tic
	bcTypes = [bcTypes_interior, bcTypes_exterior];
	dom = GalerkinAssembler2d_rxndiff.assembleBoundary(dom_geo,bcTypes,BCs); 
	dom = GalerkinAssembler2d_rxndiff.assembleMesh(dom,p,base);
	dom = GalerkinAssembler2d_rxndiff.assembleTimeStepping(dom,T,dt);
	auxfun = GalerkinAssembler2d_rxndiff.assembleCoefficients(c,k,r,f,u_o);
executionTime = toc;
fprintf(' %f s\n',executionTime)

% Solve
fprintf(' Solving:'), tic
	%prob = NewtonGalerkinSolver2d_rxndiff(dom,auxfun);
	prob2 = GalerkinSolver2d_heat(dom,auxfun);
executionTime = toc; 
fprintf(' %f s\n',executionTime)

