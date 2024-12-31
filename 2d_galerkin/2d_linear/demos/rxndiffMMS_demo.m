% rxndiffMMS_demo
clear all; x = sym('x',[1 2],'real'); syms t;
% USER INPUTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% domain bounds
%xLim_dom = [0 3];
xLim_dom = [0 1];
yLim_dom = [0 1];

% Y bounds
xLim_Y = [0 1];
yLim_Y = [0 1];

% time stepping
T = 1;

% number of inclusions
eps = 1/2;
incRatio = pi/2; % <~~~ incRatio = |delta Q| / |Y|

% mms parameters
base = 2;
demo = 0;

% specify BCs
bTypes_outer = 'TTTT';
%bTypes_outer = 'RRRR';
bTypes_inner = 'R';

% specify coefficients
p = 1 + x(1) * x(2) * t;
k = 1 + x(1) * x(2) * t;
r = 1 + x(1) * x(2) * t; 

% specify desired result
uTrue = sin(pi/2 * x(1)) * sin(pi/2 * x(2)) * t + t;


% MMS TEST %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('MMS Test Begun\n')

% assemble inputs
auxfun    = ManufacturedFunctions2d_rxndiff(p,k,r,uTrue);
mmsparams = MMSParams(base,demo=demo,timeOffset=3,timeFactor=2,pmin=3,pmax=5, ...
				meshInclusions='on',effectiveRegion='Omega_eps');

% build domain
fprintf('Initialization\n')
fprintf(' Contructing Domain:'), tic
	%dom = Domain2d(xLim_dom,yLim_dom);
	%inc = Inclusion2d_circle(xLim_Y,yLim_Y,incRatio);
	inc = Inclusion2d_square(xLim_Y,yLim_Y,incRatio);
	dom = Domain2d_punctured(xLim_dom,yLim_dom,inc,eps);
	dom = dom.add_yline(0.5);
	dom = dom.setBCTypes([bTypes_outer,bTypes_inner]);
	dom.time = TimeStepping(T,1);
executionTime = toc; 
fprintf(' %f s\n',executionTime)


% run mms test
if demo == 0
	mms = GalerkinMMS2d_rxndiff(dom,auxfun,mmsparams,errType="Linfty(L2)")

% run demo test
else
	mms = GalerkinMMS2d_rxndiff(dom,auxfun,mmsparams);
	prob = mms.problems{1};
end
