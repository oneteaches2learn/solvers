% DEMO2_DRIVER 
clear all; x = sym('x',[1 2],'real'); syms t;
% USER INPUTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% domain bounds
xLim = [0 1];
yLim = [0 1];

% time stepping
T = 1;

% number of inclusions
N_x = 1;
N_y = 1;
alpha = 2; % <~~~ alpha = |delta Q| / |Y|

% generate mesh
p = 1;
base = 2;
h = base^-p;

% specify BCs
bTypes = {'R' 'R' 'R' 'R'};
bTypes2 = 'R';

% specify coefficients
p     = 1;
k     = 1;
r     = 1;
uStar = 1;

% specify desired result
uTrue = sin(pi/2 * x(1)) * sin(pi/2 * x(2)) * t + t;

mms_test = 1;


% MMS TEST %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('MMS Test Begun\n')

% assemble inputs
bound  = PuncturedBoundary2d(bTypes,{@()(0.0),@()(0.0),@()(0.0),@()(0.0)},bTypes2,{@()(0.0)});
auxfun = PennesAuxFunctions2d(p,k,r,uStar,@()(0.0),uTrue);
time   = TimeStepping(T,1);

% build domain
fprintf('Initialization\n')
fprintf(' Contructing Domain:'), tic
incMod = InclusionModule1(alpha);
dom    = PuncturedDomain2d(xLim,yLim,N_x,N_y,incMod);
%dom    = Domain2d(xLim,yLim);
dom = dom.setEdgeBCTypes(bound);
executionTime = toc; 
fprintf(' %f s\n',executionTime)

% run mms test
if mms_test == 1
	mms = GalerkinPennes2d_mms(dom,time,auxfun,errType="Linfty(L2)")

% run demo test
elseif mms_test == 0
	mms = GalerkinPennes2d_mms(dom,time,auxfun,demo=p);
	prob = mms.problems{1};
end
