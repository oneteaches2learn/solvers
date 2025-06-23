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
%hand_small
hand_small_fine

% time stepping
T  = 6;
dt = 0.025;

% PDE INFORMATION
% specify physical conditions
u_body = 37;
u_air = -40;

% specify coefficients
c = 10^3;
k = 0.5;
f = 0;

% specify boundary conditions
bcTypes_exterior = 'RRRR';
bcTypes_interior = 'R';
bc_air = {1,u_air};
bc_wrist = {100,v};
BCs = {bc_wrist,bc_air,bc_air,bc_air,{0,0}};

% specify initial condition
u_o = u_body;


% REACTION TERM
% specify reaction term
r_const = 10^2;

%{
% linear model
r_activ = 1;
r_activ_du = 0;
%}

% ramp activation model
v_min = 28; v_max = 37; u_min = 10; u_max = 37; gamma = 0;
r_activ = @(x1,x2,t,u,v) Coefficients.ramp_activation(u,v, ...
							u_min=u_min, u_max=u_max, v_min=v_min, v_max=v_max, gamma=gamma);
r_activ_du = @(x1,x2,t,u,v) Coefficients.ramp_activation_du(u,v, ...
							u_min=u_min, u_max=u_max, v_min=v_min, v_max=v_max, gamma=gamma);

%{
% logistic activation model
u_L = 2; u_k = 0.2; u_0 = 32;
v_L = 3; v_k = 0.5; v_0 = 38.4;
r_activ = @(x1,x2,t,u,v) Coefficients.logistic_activation(u,v, ...
							u_L=u_L, u_k=u_k, u_0=u_0, v_L=v_L, v_k=v_k, v_0=v_0);
r_activ_du = @(x1,x2,t,u,v) Coefficients.logistic_activation_du(u,v, ...
							u_L=u_L, u_k=u_k, u_0=u_0, v_L=v_L, v_k=v_k, v_0=v_0);
%}

% ODE INFORMATION
g = 1;
v_o = u_body;
s = 1; % <~~~ NOTE: This is a placeholder value
s_const = 10^0;
constraints.vLower = NaN;
constraints.vUpper = 37;
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
	dom = GalerkinAssembler2d_frostbite.assembleBoundary(dom_geo,bcTypes,BCs); 
	dom = GalerkinAssembler2d_frostbite.assembleNonlinearBoundary(dom,bcTypes,BCs);
	dom.boundary = dom.boundary.setBoundaryNodeLists(dom.mesh);
	dom.boundary = dom.boundary.setFreeNodes(dom.mesh);
	dom = GalerkinAssembler2d_frostbite.assembleTimeStepping(dom,T,dt);
	auxfun = GalerkinAssembler2d_frostbite.assembleCoefficients(c,k,f,u_o);
	auxfun = GalerkinAssembler2d_frostbite.assembleReactionTerm(auxfun,r_const,r_activ,r_activ_du);
executionTime = toc;
fprintf(' %f s\n',executionTime)

% Assemble ODE Data
fprintf(' Assembling ODE Data:'), tic
	data.g = g;
	data.vInit = v_o;
	data.cofs.s = s;
	data.time.T = T;
	data.time.dt = dt;
	data.constraints = constraints;
	options.order = order;
	ode = ODE(data,options);
executionTime = toc;
fprintf(' %f s\n',executionTime)

% Set s for ODE
r_activ = @(u,v) auxfun.cofs.r_activ(0,0,0,u,v);
data.cofs.s = @(u,v) s_const * r_activ(u,v);


% Solve
fprintf(' Solving:'), tic
	ode = ODE(data,options);
	%prob = CoupledNewtonSolver2d_rxndiff(dom,auxfun,ode);
	prob = CoupledNewtonSolver2d_frostbite(dom,auxfun,ode);
executionTime = toc; 
fprintf(' %f s\n',executionTime)
