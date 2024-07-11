% DEMO2_DRIVER 
clear all; x = sym('x',[1 2],'real');
% USER INPUTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% domain bounds
xLim = [0 1];
yLim = [0 1];

% number of inclusions
N_x = 1;
N_y = 1;
alpha = 3;

% generate mesh
p = 3;
base = 2;
h = base^-p;

% specify BCs
bTypes = {'R' 'R' 'R' 'R'};
bTypes2 = 'R';

% specify coefficients
k = 1;
r = 1;

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
incMod = InclusionModule1(alpha);
dom    = PuncturedDomain2d(xLim,yLim,N_x,N_y,incMod);

dom.plotGeometry


