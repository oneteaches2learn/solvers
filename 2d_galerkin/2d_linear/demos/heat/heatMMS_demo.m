% HeatMMS_demo
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
incRatio = pi/2; % <~~~ incRatio = |delta Q| / |Y|

% mms parameters
base = 2;
demo = 0;

% specify BCs
bTypes_outer = 'RRRR';
bTypes_inner = 'N';

% specify coefficients
%c = 2 + sin(x(1) + x(2)) * t + t^2; 
%k = 2 + sin(x(1) + x(2)) * t + t^2; 
c = 1;
k = 1;


% specify desired result
%uTrue = sin(2 * pi * x(1)) * sin(2 * pi * x(2)) * (t + 1) + (t + 1);
uTrue = cos(2 * pi * x(1)) * cos(2* pi * x(2));


% MMS TEST %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('MMS Test Begun\n')

% assemble inputs
auxfun    = ManufacturedFunctions2d_heat(c,k,uTrue);
mmsparams = MMSParams(base,demo=demo,timeOffset=2,timeFactor=2,pmin=2,pmax=4, ...
				meshInclusions='on',effectiveRegion='Omega_eps');

% build domain
fprintf('Initialization\n')
fprintf(' Contructing Domain:'), tic
	dom = Domain2d(xLim_dom,yLim_dom);
	%inc = Inclusion2d_circle(xLim_Y,yLim_Y,incRatio);
	%inc = Inclusion2d_square(xLim_Y,yLim_Y,incRatio);
	%dom = Domain2d_punctured(xLim_dom,yLim_dom,inc,eps);
	%dom = dom.add_yline(0.5);
	dom = dom.setBCTypes([bTypes_outer,bTypes_inner]);
	dom.time = TimeStepping(T,1);
executionTime = toc; 
fprintf(' %f s\n',executionTime)

% run mms test
if demo == 0
	mms = GalerkinMMS2d_heat(dom,auxfun,mmsparams,errType="Linfty(L2)")

% run demo test
else
	mms = GalerkinMMS2d_heat(dom,auxfun,mmsparams);
	prob = mms.problems{1};
end
