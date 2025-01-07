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
T  = 6;
dt = 0.05;

% PDE INFORMATION
% specify physical conditions
u_body = 37;
u_air = -40;

% specify coefficients
c = 10^3;
k = 0.5;
r = 1; % <~~~ NOTE: This is a placeholder value
r_const = 10^2;
f = 0;

% specify boundary conditions
bcTypes_exterior = 'RRRR';
bcTypes_interior = 'R';
bc_air = {1,u_air};
bc_wrist = {100,v};
BCs = {bc_wrist,bc_air,bc_air,bc_air,{0,0}};

% specify initial condition
u_o = u_body;


% ODE INFORMATION
g = 1;
v_o = u_body;
s = 1; % <~~~ NOTE: This is a placeholder value
s_const = 10^-1;
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


% TEMPORARY: Manually edit coefficients
v_min = 28;
v_max = 37;
u_min = 10;
u_max = 32;
r_cof = @(u,v) Coefficients.physiological_coefficient1(u,v, ...
							u_min=u_min, u_max=u_max, v_min=v_min, v_max=v_max);
r_couple = @(u,v) (u - v);
dr_du = @(u,v) Coefficients.physiological_coefficient1_du(u,v, ...
							u_min=u_min, u_max=u_max, v_min=v_min, v_max=v_max);
auxfun.cofs.r = @(x1,x2,u,v) r_const * r_cof(u,v) .* r_couple(u,v);
auxfun.cofs.dr_du = @(x1,x2,u,v) r_const * dr_du(u,v);

data.cofs.s = @(u,v) s_const * r_cof(u,v);
ode = ODE(data,options);

% Solve
fprintf(' Solving:'), tic
	prob = CoupledNewtonSolver2d_rxndiff(dom,auxfun,ode);
executionTime = toc; 
fprintf(' %f s\n',executionTime)
