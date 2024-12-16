% rxndiffSolve_demo 
clear all; x = sym('x',[1 2],'real'); syms t; syms u;
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

% domain parameters
base = 2;
p_start = 3;
p_end = 5;
T = 1;

% specify coefficients
c = 1;
k = 1;
r = u^2 - u;

% specify boundary conditions
bcTypes_interior = 'DDDD';
bcTypes_exterior = 'D';

% desired solution
uTrue = sin(pi/2 * x(1)) * sin(pi/2 * x(2)) * t + t^2;

% BLACK BOX %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SETUP TRIAL
fprintf('RxnDiff MMS Test Begun\n')

fprintf(' Contructing Domain Geometry:'), tic
	%inc = Inclusion2d_square(xLim_Y,yLim_Y,incRatio);
	inc = Inclusion2d_circle(xLim_Y,yLim_Y,incRatio);
	dom_geo = GalerkinAssembler2d_rxndiff.assembleDomainGeometry(xLim_dom,yLim_dom,inc,eps);
	%dom_geo = GalerkinAssembler2d_rxndiff.assembleDomainGeometry(xLim_dom,yLim_dom);
executionTime = toc; 
fprintf(' %f s\n',executionTime)

fprintf(' Manufacturing Data:'), tic
	bcTypes = [bcTypes_interior, bcTypes_exterior];
	manfun = ManufacturedFunctions2d_rxndiff(c,k,r,uTrue);
	manfun = manfun.manufactureBCs(dom_geo,bcTypes);
executionTime = toc; 
fprintf(' %f s\n',executionTime)

% collect coefficients
c = manfun.c;
k = manfun.k;
r = manfun.r;
f = manfun.f;
u_o = manfun.uInit;
BCs = manfun.bcConds;


% RUN MMS
for p = p_start:p_end

	i = (p - p_start) + 1;
	dt = T / base^(2*(p-p_start));

	fprintf(' Problem 1:\n')
	fprintf('  Assembly:'), tic
		dom = GalerkinAssembler2d_rxndiff.assembleBoundary(dom_geo,bcTypes,BCs); 
		dom = GalerkinAssembler2d_rxndiff.assembleMesh(dom,p,base);
		dom = GalerkinAssembler2d_rxndiff.assembleTimeStepping(dom,T,dt);
		auxfun = GalerkinAssembler2d_rxndiff.assembleCoefficients(c,k,r,f,u_o);
	executionTime = toc;
	fprintf(' %f s\n',executionTime)

	fprintf('  Solving:'), tic
		prob{i} = NewtonGalerkinSolver2d_rxndiff(dom,auxfun);
	executionTime = toc; 
	fprintf(' %f s\n',executionTime)

end


% ANALYZE RESULTS
% store variables
trials = length(prob);

% Store errors
errors = zeros(1,trials);
ratios = zeros(1,trials-1);
orders = zeros(1,trials-1);

% Compute error
fprintf('Computing Errors:'), tic
for i = 1:trials

	% compute quadrature (on each timestep, if time-varying)
	sol = prob{i}.solution;
	err = prob{i}.domain.L2err_centroidQuadrature_nodal(sol,uTrue);

	% compute L_infty error across all time steps
	errors(i) = max(err);

end

% Compute ratios and orders
for i = 2:trials
	ratios(i-1) = errors(i-1)/errors(i);
	orders(i-1) = log(ratios(i-1)) / log(base);
end
executionTime = toc;
fprintf('%f s\n',executionTime) 

orders