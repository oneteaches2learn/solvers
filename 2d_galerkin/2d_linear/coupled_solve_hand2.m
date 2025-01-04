% rxndiffSolve_demo 
clear all; x = sym('x',[1 2],'real'); syms t; syms u; syms v;
% USER INPUTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% DOMAIN INFORMATION
% domain shape
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

% time stepping
T  = 3;
dt = 0.05;

% PDE INFORMATION
% specify physical conditions
u_body = 37;
u_air = -40;

% specify coefficients
c = 10^3;
k = 0.5;
r = (u - v);
f = 0;

% specify boundary conditions
bcTypes_exterior = 'RRRR';
bcTypes_interior = 'R';
bc_air = {1,u_air};
bc_wrist = {1,v};
BCs = {bc_wrist,bc_air,bc_air,bc_air,{0,0}};

% specify initial condition
u_o = u_body;


% ODE INFORMATION
g = 0;
v_o = u_body;
s = 1;
order = 1;


% BLACK BOX %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SETUP TRIAL
fprintf('Solving Coupled System\n')

% Construct Domain Geometry
fprintf(' Contructing Domain Geometry:'), tic
	dom_geo = Domain2d.domainFromGmsh(msh);
executionTime = toc; 
fprintf(' %f s\n',executionTime)


% Assemble PDE Data 
fprintf(' Assembling PDE Data:'), tic
	bcTypes = [bcTypes_exterior, bcTypes_interior];
	dom = GalerkinAssembler2d_rxndiff.assembleBoundary(dom_geo,bcTypes,BCs); 
	dom = GalerkinAssembler2d_rxndiff.assembleNonlinearBoundary(dom,bcTypes,BCs);
	dom.boundary = dom.boundary.setBoundaryNodeLists(dom.mesh);
	dom.boundary = dom.boundary.setFreeNodes(dom.mesh);
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
	prob = CoupledNewtonSolver2d_rxndiff(dom,auxfun,ode);
executionTime = toc; 
fprintf(' %f s\n',executionTime)


% AUXILIARY FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y = ramp_function(x)

	y = x;
	y(x<0) = 0;

end