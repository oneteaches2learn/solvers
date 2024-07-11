% DEMO2_DRIVER 
clear all; x = sym('x',[1 2],'real');
% USER INPUTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% domain bounds
xLim = [0 1];
yLim = [0 1];

% number of inclusions
N_x = 3;
N_y = 3;
alpha = 1; % <~~~ alpha = |delta Q| / |Y|

% generate mesh
p = 1;
base = 2;
h = base^-p;

% specify BCs
bTypes = {'R' 'R' 'R' 'R'};
bTypes2 = 'R';

% specify coefficients
k = 3 + x(1)*x(2);
r = 1 + x(1) * x(2);

% specify desired result
uTrue = sin(pi/2*x(1))*sin(pi/2*x(2));

mms_test = 1;


% MMS TEST %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('MMS Test Begun\n')
% ensure functions have correct variables
uTrue = symfun(uTrue,x);

% assemble inputs
bound  = PuncturedBoundary2d(bTypes,{@()(0.0),@()(0.0),@()(0.0),@()(0.0)},bTypes2,{@()(0.0)});
auxfun = PoissonAuxFunctions2d(k,r,uTrue);

% build domain
fprintf('Initialization\n')
fprintf(' Contructing Domain:'), tic
incMod = InclusionModule1(alpha);
dom    = PuncturedDomain2d(xLim,yLim,N_x,N_y,incMod);
%dom    = Domain2d(xLim,yLim);
dom = dom.setEdgeBCTypes(bound);
executionTime = toc; 
fprintf(' %f s\n',executionTime)


if 1 == 1
% run mms test
if mms_test == 1
	mms = GalerkinPoisson2d_mms(dom,auxfun,errType="L2")

% run demo test
elseif mms_test == 0
	mms = GalerkinPoisson2d_mms(dom,auxfun,demo=p);
	prob = mms.problems{1};
end
end
