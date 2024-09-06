% volume_fraction
% Test to compute the differential equation that gives the coefficients of the
% tensor for the stiffness matrix
clear all; x = sym('x',[1 2],'real'); syms t;
% USER INPUTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% domain bounds
xLim_dom = [0 1];
yLim_dom = [0 1];

% Y bounds
xLim_Y = [0 1];
yLim_Y = [0 1];

% number of inclusions
eps = 1;
incRatio = pi/2; % <~~~ incRatio = |delta Q| / |Y|

% mesh parameters
p = 5;
base = 2;

% specify source
f = 0;

% specify BCs
bTypes = {'P' 'P' 'P' 'P'};
bTypes2 = 'N';

% specify coefficients
k = 1;
r = 0;

% specify neumann BC
dy1_dn = {0,1,0,-1};
dy2_dn = {-1,0,1,0};


% MMS TEST %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Test Begun\n')

% assemble inputs
bound     = Boundary2d_punctured(bTypes,{@()(0.0),@()(0.0),@()(0.0),@()(0.0)},bTypes2,dy1_dn);

% build domain
fprintf(' Contructing Domain:'), tic
	%inc = Inclusion2d_circle(xLim_Y,yLim_Y,incRatio);
	inc = Inclusion2d_square(xLim_Y,yLim_Y,incRatio);
	dom = GalerkinPoisson2d_assembler.assembleDomainGeometry(xLim_dom,yLim_dom,inc,eps);
	dom = GalerkinPoisson2d_assembler.assembleBoundary(dom,bTypes,0,u_N,0,0,bTypes2); 
	dom = GalerkinPoisson2d_assembler.assembleMesh(dom,p,base); 
executionTime = toc; 
fprintf(' %f s\n',executionTime)

% assemble other parameters
fprintf(' Assembling Parameters:'), tic
	cofs   = GalerkinPoisson2d_assembler.assembleCoefficients(c,k);
	uInit  = GalerkinPoisson2d_assembler.assembleInitialCondition(u_o);
	time   = GalerkinPoisson2d_assembler.assembleTimeStepping(T,dt);
	source = GalerkinPoisson2d_assembler.assembleSource(f);
executionTime = toc; 
fprintf(' %f s\n',executionTime)

% RUN TRIAL
fprintf(' Generating Solution:'), tic
	prob = GalerkinPoisosn2d_solver(dom,time,cofs,uInit,source);
executionTime = toc; 
fprintf(' %f s\n',executionTime)