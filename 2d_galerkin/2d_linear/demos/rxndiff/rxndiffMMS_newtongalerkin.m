% PoissonMMS_demo
clear all; x = sym('x',[1 2],'real'); syms t; syms u;
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
Eps = 1/2;
incRatio = pi / 2; % <~~~ incRatio = |delta Q| / |Y|
effRegion = 'Omega_eps';
meshInclusions = 'off';

% mms parameters
base = 2;
demo = 0;

% specify BCs
bTypes_outer = 'NNNN';
bTypes_inner = 'R';

% specify coefficients
%{
c = 1 + x(1) + x(2);
k = 1 + x(1) + x(2) + t;
r = (u^4 - 1);
r = 0;
%}
c = 1;
k = 1;
r = u^3;

% specify nonlinear boundary conditions
%u_N = sin(u); 
u_N = 1*(u + u^3);
%u_R = u^3;
%u_R = u^2;
%alpha_R = 2 + x(1) * x(2);
%u_R = exp(u);
%alpha_R = 1 + x(1) + x(2) + t;
u_R = u^3;
alpha_R = 1;

% specify desired result
%uTrue = 1 + sin(pi / 2* x(1)) * sin(pi / 2 * x(2));
%uTrue = cos(2 * pi * x(1)) * cos(2* pi * x(2));
%uTrue = exp(x(1) * x(2));
%uTrue = sin(pi * x(1)) * sin(pi * x(2)) * t;
uTrue = sin(pi / 2 * x(1)) * sin(pi / 2 * x(2)) * t + t + 1;
%uTrue = 1 + t;
%uTrue = 1;


% MMS TEST %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('MMS Test Begun\n')

% assemble inputs
auxfun    = ManufacturedFunctions2d_rxndiff(c,k,r,uTrue,u_N=u_N,alpha_R=alpha_R,u_R=u_R);
mmsparams = MMSParams(base,demo=demo,timeOffset=2,timeFactor=2,pmin=4,pmax=6, ...
				meshInclusions=meshInclusions,effectiveRegion=effRegion);

% build dom_eps_epsain
fprintf('Initialization\n')
fprintf(' Contructing Domain:'), tic
	dom = Domain2d(xLim_dom,yLim_dom);
	%inc = Inclusion2d_circle(xLim_Y,yLim_Y,incRatio);
	%inc = Inclusion2d_square(xLim_Y,yLim_Y,incRatio);
	%dom = Domain2d_punctured(xLim_dom,yLim_dom,inc,Eps);
	%dom = dom.add_yline;
	dom = dom.setBCTypes([bTypes_outer,bTypes_inner]);
	dom.time = TimeStepping(T,1);
executionTime = toc; 
fprintf(' %f s\n',executionTime)

% run mms test
if demo == 0
	mms = NewtonGalerkinMMS2d_rxndiff(dom,auxfun,mmsparams,errType="L2")

% run demo test
else
	mms = NewtonGalerkinMMS2d_rxndiff(dom,auxfun,mmsparams);
		prob = mms.problems{1};
end