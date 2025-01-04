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
T  = 5;
dt = 0.01;

% specify physical conditions
u_body = 37;
u_air = -40;

% specify coefficients
c = 10^3;
k = 0.5;
r = 0;
uStar = u_body;
f = 0;

% specify boundary conditions
bcTypes_interior = 'RRRR';
bcTypes_exterior = 'R';
D = 1;
R_body = {100,u_body};
R_air = {5,u_air};
BCs = {R_body,R_air,R_air,R_air,D};

% specify initial condition
u_o = u_body;


% BLACK BOX %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SETUP TRIAL
fprintf('Solving Nonlinear RxnDiff\n')

% Construct Domain Geometry
fprintf(' Contructing Domain Geometry:'), tic
	%domain_dist50
	%rect_coarse
	%fist2
	%diamond_backwards	
	%rect
	%rect2
	%diamond2
	%skew
	hand_small
	%hand_small_fine
	dom_geo = Domain2d.domainFromGmsh(msh);
executionTime = toc; 
fprintf(' %f s\n',executionTime)


% Assemble Data
fprintf(' Assembling Data:'), tic
	bcTypes = [bcTypes_interior, bcTypes_exterior];
	%dom = dom_geo;
	dom = GalerkinAssembler2d_rxndiff.assembleBoundary(dom_geo,bcTypes,BCs); 
	dom.boundary = dom.boundary.setBoundaryNodeLists(dom.mesh);
	dom.boundary = dom.boundary.setFreeNodes(dom.mesh);
	dom = GalerkinAssembler2d_rxndiff.assembleTimeStepping(dom,T,dt);
	auxfun = GalerkinAssembler2d_rxndiff.assembleCoefficients(c,k,r,f,u_o);
	auxfun = GalerkinAssembler2d_pennes.assembleCoefficients(c,k,r,uStar,f,u_o);
	%auxfun.cofs.uStar = @(x,y) uStar;
executionTime = toc;
fprintf(' %f s\n',executionTime)

% Solve
fprintf(' Solving:'), tic
	%prob = NewtonGalerkinSolver2d_rxndiff(dom,auxfun);
	prob = GalerkinSolver2d_pennes(dom,auxfun);
executionTime = toc; 
fprintf(' %f s\n',executionTime)

