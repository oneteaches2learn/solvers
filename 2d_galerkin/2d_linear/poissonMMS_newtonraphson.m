% PoissonMMS_demo
clear all; x = sym('x',[1 2],'real'); syms t; syms u;
% USER INPUTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% domain bounds
xLim_dom = [0 1];
yLim_dom = [0 1];

% Y bounds
xLim_Y = [0 1];
yLim_Y = [0 1];

% number of inclusions
Eps = 1/2;
incRatio = 1; % <~~~ incRatio = |delta Q| / |Y|
effRegion = 'Omega_eps';
meshInclusions = 'off';

% mms parameters
base = 2;
demo = 0;

% specify BCs
bTypes_outer = 'RRRR';
bTypes_inner = 'R';

% specify coefficients
k = 1;
r = u^3 - u;

% specify desired result
uTrue = sin(pi * x(1)) * sin(pi * x(2));


% MMS TEST %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('MMS Test Begun\n')

% assemble inputs
auxfun    = ManufacturedFunctions2d_poisson(k,r,uTrue);
mmsparams = MMSParams(base,demo=demo,timeOffset=4,timeFactor=2,pmin=4,pmax=6, ...
				meshInclusions=meshInclusions,effectiveRegion=effRegion);

% build dom_eps_epsain
fprintf('Initialization\n')
fprintf(' Contructing Domain:'), tic
	%dom = Domain2d(xLim_dom,yLim_dom);
	inc = Inclusion2d_circle(xLim_Y,yLim_Y,incRatio);
	%inc = Inclusion2d_square(xLim_Y,yLim_Y,incRatio);
	dom = Domain2d_punctured(xLim_dom,yLim_dom,inc,Eps);
	dom = dom.add_yline;
	dom = dom.setBCTypes([bTypes_outer,bTypes_inner]);
executionTime = toc; 
fprintf(' %f s\n',executionTime)

% run mms test
if demo == 0
	mms = NewtonRaphsonMMS2d_poisson(dom,auxfun,mmsparams,errType="L2")

% run demo test
else
	mms = NewtonRaphsonMMS2d_poisson(dom,auxfun,mmsparams);
	prob = mms.problems{1};
end