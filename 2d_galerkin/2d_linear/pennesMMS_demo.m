% PENNESMMS_DEMO
clear all; x = sym('x',[1 2],'real'); syms t;
% USER INPUTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% domain bounds
xLim_dom = [0 2];
yLim_dom = [0 1];

% Y bounds
xLim_Y = [0 1];
yLim_Y = [0 1];

% time stepping
T = 1;

% number of inclusions
eps = 1;
incRatio = pi/2; % <~~~ incRatio = |delta Q| / |Y|

% mms parameters
base = 2;
demo = 0;

% specify BCs
bTypes = {'R' 'N' 'D' 'D'};
bTypes2 = 'R';

% specify coefficients
p = 1 + x(1) * x(2) * t;
k = 1 + x(1) * x(2) * t;
r = 1 + x(1) * x(2) + t; 
uStar = 1 + x(1) * x(2) + t;

% specify desired result
uTrue = sin(pi/2 * x(1)) * sin(pi/2 * x(2)) * t;


% MMS TEST %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('MMS Test Begun\n')

% assemble inputs
bound     = Boundary2d_punctured(bTypes,{@()(0.0),@()(0.0),@()(0.0),@()(0.0)},bTypes2,{@()(0.0)});
auxfun    = ManufacturedFunctions2d_pennes(p,k,r,uStar,uTrue);
mmsparams = MMSParams(base,demo=demo,timeOffset=3,timeFactor=2,pmin=3,pmax=5);

% build domain
fprintf('Initialization\n')
fprintf(' Contructing Domain:'), tic
	inc = Inclusion2d_circle(xLim_Y,yLim_Y,incRatio);
	%inc = Inclusion2d_square(xLim_Y,yLim_Y,incRatio);
	dom = Domain2d_punctured(xLim_dom,yLim_dom,inc,eps);
	dom = dom.setEdgeBCTypes(bound);
	dom.time = TimeStepping(T,1);
executionTime = toc; 
fprintf(' %f s\n',executionTime)

% run mms test
if demo == 0
	mms = GalerkinMMS2d_pennes(dom,auxfun,mmsparams,errType="Linfty(L2)")

% run demo test
else
	mms = GalerkinMMS2d_pennes(dom,auxfun,mmsparams);
	prob = mms.problems{1};
end