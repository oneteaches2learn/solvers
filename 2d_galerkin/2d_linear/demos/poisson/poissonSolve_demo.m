% PoissonMMS_demo
clear all; x = sym('x',[1 2],'real'); syms t;
% USER INPUTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% domain bounds
xLim_dom = [0 1];
yLim_dom = [0 1];

% Y bounds
xLim_Y = [0 1];
yLim_Y = [0 1];

% number of inclusions
epsilon = 1;
incRatio = pi / 2; % <~~~ incRatio = |delta Q| / |Y|

% mms parameters
p = 2;
base = 2;

% specify coefficients
k = 2 + sin(x(1) + x(2));
r = 2 + sin(x(1) + x(2));

% specify source
f = 1;

% specify BCs
bTypes_outer = 'DDDD';
bTypes_inner = 'R';

% specify Dirichlet conditions
u_D = 0;

% specify Neumann conditions
u_N = 1;

% specify Robin conditions
beta = 1;
u_R = 1;



% MMS TEST %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Poisson Trial Begun\n')

% assemble domain
fprintf(' Contructing Domain:'), tic
	inc = Inclusion2d_circle(xLim_Y,yLim_Y,incRatio);
	%inc = Inclusion2d_square(xLim_Y,yLim_Y,incRatio);
	%dom = GalerkinAssembler2d_poisson.assembleDomainGeometry(xLim_dom,yLim_dom,inc,epsilon);
	%dom = GalerkinAssembler2d_poisson.assembleDomainGeometry(xLim_dom,yLim_dom);
	dom = Domain2d(xLim_dom,yLim_dom);
	%dom = dom.add_yline;
	%dom = dom.inclusionsOFF;
	dom = GalerkinAssembler2d_poisson.assembleBoundary(dom,bTypes_outer,u_D,u_N,beta,u_R,bTypes_inner); 
	dom = GalerkinAssembler2d_poisson.assembleMesh(dom,p,base); 
executionTime = toc; 
fprintf(' %f s\n',executionTime)

% assemble other parameters
fprintf(' Assembling Parameters:'), tic
	auxfun = GalerkinAssembler2d_poisson.assembleCoefficients(k,r,f);
executionTime = toc; 
fprintf(' %f s\n',executionTime)

% run trial
fprintf(' Generating Solution:'), tic
	prob = GalerkinSolver2d_poisson(dom,auxfun);
executionTime = toc; 
fprintf(' %f s\n',executionTime)