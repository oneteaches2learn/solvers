% packageMMS_demo
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
eps = 1;
incRatio = 1; % <~~~ incRatio = |delta Q| / |Y|

% specify BCs
bTypes = {'R' 'R' 'R' 'R'};
bTypes2 = 'R';

% specify coefficients
p = 1 + x(1) * x(2) * t;
k = 1 + x(1) * x(2) * t;
r = 1 + x(1) * x(2) * t; 
uStar = 1 + x(1) * x(2) * t;

% specify desired result
uTrue = sin(pi/2 * x(1)) * sin(pi/2 * x(2)) * t + t;


% ADDN'L INPUTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% mms parameters
base = 2;
demo = 0;
timeOffset = 1;
timeFactor = 2;


% MMS TEST %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('MMS Test Begun\n')
fprintf(' Expected Runtime: 1 to 2 minutes\n')
fprintf(' First run will take ~30s longer because parallel pool must be initiated.\n')

% assemble inputs
bound     = Boundary2d_punctured(bTypes,{@()(0.0),@()(0.0),@()(0.0),@()(0.0)},bTypes2,{@()(0.0)});
auxfun    = ManufacturedFunctions2d_pennes(p,k,r,uStar,uTrue);
time      = TimeStepping(T,1);
mmsparams = MMSParams(base,demo=demo,timeOffset=timeOffset,timeFactor=timeFactor,pmin=2,pmax=5);

% build domain
fprintf('Initialization\n')
fprintf(' Contructing Domain:'), tic
	inc = Inclusion2d_circle(xLim_Y,yLim_Y,incRatio);
	%inc = Inclusion2d_square(xLim_Y,yLim_Y,incRatio);
	dom = Domain2d_punctured(xLim_dom,yLim_dom,inc,eps);
	dom = dom.setEdgeBCTypes(bound);
executionTime = toc; 
fprintf(' %f s\n',executionTime)


% run mms test
if demo == 0
	mms = GalerkinPennes2d_mms(dom,time,auxfun,mmsparams,errType="Linfty(L2)")

% run demo test
else
	mms = GalerkinPennes2d_mms(dom,time,auxfun,mmsparams);
	prob = mms.problems{1};
end
