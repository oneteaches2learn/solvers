clear all; x = sym('x',[1 2],'real'); syms t;

% Define params structure
params.xLim_dom = [0, 1];
params.yLim_dom = [0, 1];
params.xLim_Y = [0, 1];
params.yLim_Y = [0, 1];
params.T = 1;
params.incRatio = 1;
params.eps = 1/2;
params.bTypes_outer = 'NDRR';
params.bTypes_inner = 'R';
params.k = 2 + sin(x(1) + x(2));
params.r = 2 + sin(x(1) + x(2));
params.demo = 0;
params.uTrue = cos(2 * pi * x(1)) * cos(2 * pi * x(2));

options.base = 2;
options.demo = 0;
options.timeOffset = 4;
options.timeFactor = 2;
options.pmin = 4;
options.pmax = 8;
options.meshInclusions = 'on';
options.effectiveRegion = 'Omega_eps';


% Call the function
poissonMMS_test(params,options)