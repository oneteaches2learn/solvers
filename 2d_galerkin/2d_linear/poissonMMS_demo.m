% PoissonMMS_demo
clear all; x = sym('x',[1 2],'real'); syms t;
% USER INPUTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% domain bounds
xLim_dom = [0 1];
yLim_dom = [0 1];

% Y bounds
xLim_Y = [0 1];
yLim_Y = [0 1];

% time stepping
T = 1;

% number of inclusions
eps = 1/2;
incRatio = 1; % <~~~ incRatio = |delta Q| / |Y|

% mms parameters
base = 2;
demo = 0;

% specify BCs
bTypes_outer = 'NDRR';
bTypes_inner = 'R';

% specify coefficients
k = 2 + sin(x(1) + x(2));
r = 2 + sin(x(1) + x(2));

% specify desired result
%uTrue = sin(pi/2 * x(1)) * sin(pi/2 * x(2));
%uTrue = sin(2 * pi * x(1)) * sin(2 * pi * x(2));
uTrue = cos(2 * pi * x(1)) * cos(2 * pi * x(2));


% MMS TEST %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('MMS Test Begun\n')

% assemble inputs
auxfun    = ManufacturedFunctions2d_poisson(k,r,uTrue);
mmsparams = MMSParams(base,demo=demo,timeOffset=4,timeFactor=2,pmin=4,pmax=8, ...
				meshInclusions='on',effectiveRegion='Omega_eps');

% build dom_eps_epsain
fprintf('Initialization\n')
fprintf(' Contructing Domain:'), tic
	%dom = Domain2d(xLim_dom,yLim_dom);
	inc = Inclusion2d_circle(xLim_Y,yLim_Y,incRatio);
	%inc = Inclusion2d_square(xLim_Y,yLim_Y,incRatio);
	dom = Domain2d_punctured(xLim_dom,yLim_dom,inc,eps);
	dom = dom.add_yline(0.8);
	dom = dom.setBCTypes([bTypes_outer,bTypes_inner]);
executionTime = toc; 
fprintf(' %f s\n',executionTime)

% run mms test
if demo == 0
	mms = GalerkinMMS2d_poisson(dom,auxfun,mmsparams,errType="L2")

% run demo test
else
	mms = GalerkinMMS2d_poisson(dom,auxfun,mmsparams);
	prob = mms.problems{1};
end