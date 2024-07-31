% HEAT_DEMO
clear all; x = sym('x',[1 2],'real'); syms t;
% USER INPUTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% domain bounds
xLim = [0 1];
yLim = [0 1];

% time stepping
T = 1;

% number of inclusions
N_x = 2;
N_y = 2;
alpha = 2; % <~~~ alpha = |delta Q| / |Y|

% mms parameters
base = 2;
demo = 0;
timeOffset = 1;
timeFactor = 2;

% specify BCs
bTypes = {'R' 'R' 'R' 'R'};
bTypes2 = 'R';

% specify coefficients
p = 1 + t;
k = 1 + t;

% specify desired result
uTrue = sin(pi/2 * x(1)) * sin(pi/2 * x(2)) * t + t;



% MMS TEST %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('MMS Test Begun\n')

% assemble inputs
bound     = PuncturedBoundary2d(bTypes,{@()(0.0),@()(0.0),@()(0.0),@()(0.0)},bTypes2,{@()(0.0)});
auxfun    = ManufacturedFunctions2d_heat(p,k,uTrue);
time      = TimeStepping(T,1);
mmsparams = MMSParams(base,demo=demo,timeOffset=timeOffset,timeFactor=timeFactor,pmax=3);

% build domain
fprintf('Initialization\n')
fprintf(' Contructing Domain:'), tic
	incMod = InclusionModule1(alpha);
	dom    = PuncturedDomain2d(xLim,yLim,N_x,N_y,incMod);
	dom    = dom.setEdgeBCTypes(bound);
executionTime = toc; 
fprintf(' %f s\n',executionTime)

% run mms test
if demo == 0
	mms = GalerkinHeat2d_mms(dom,time,auxfun,mmsparams,errType="Linfty(L2)")

% run demo test
else
	mms = GalerkinHeat2d_mms(dom,time,auxfun,mmsparams);
	prob = mms.problems{1};
end
