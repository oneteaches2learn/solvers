% poissonMMS_test_BCs
% Cycles through various coefficients
clear all; x = sym('x',[1 2],'real'); syms t; start = tic;

% FIXED PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% domain bounds
params.xLim_dom = [0, 1];
params.yLim_dom = [0, 1];

% Y bounds
params.xLim_Y = [0, 1];
params.yLim_Y = [0, 1];

% number of inclusions
params.incRatio = 1;
params.eps = 1/2;

% time stepping
params.T = 1;

% boundary conditions
params.bTypes_outer = 'RRRR';
params.bTypes_inner = 'R';

% coefficients
% params.k = ...
% params.r = ...
% ... set in loop for this test

% desired true solution
params.uTrue = cos(2 * pi * x(1)) * cos(2 * pi * x(2));

% mms parameters
options.base = 2;
options.demo = 0;
options.timeOffset = 4;
options.timeFactor = 2;
options.pmin = 4;
options.pmax = 7;
options.meshInclusions = 'on';
options.effectiveRegion = 'Omega_eps';
options.expectedRate = 2;
options.convergenceTolerance = 0.1;


% EXPIREMENTAL PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% coefficients
% test 1: spatially varying
k(1) = 2 + sin(x(1) + x(2));
r(1) = 4 + 2 * sin(x(1) + x(2));

% test 2: constant
k(2) = 2;
r(2) = 4;

% test 3: zero r
k(3) = 2 + sin(x(1) + x(2));
r(3) = 0;


% RUN TEST %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% generate test parameters
testNum = 1;
loopsize(1) = length(k);
nTests = prod(loopsize);

% create domain
%dom = Domain2d(xLim_dom,yLim_dom);
inc = Inclusion2d_circle(params.xLim_Y,params.yLim_Y,params.incRatio);
%inc = Inclusion2d_square(params.xLim_Y,params.yLim_Y,params.incRatio);
dom = Domain2d_punctured(params.xLim_dom,params.yLim_dom,inc,params.eps);
dom = dom.add_yline(0.8);
dom = dom.setBCTypes([params.bTypes_outer,params.bTypes_inner]);

for i = 1:loopsize(1), 

    % store experimental variable in params
    params.k = k(i);
    params.r = r(i);

    % run test
    fprintf('TEST %d of %d. poissonMMS. coefficients: k = %s, r = %s\n', ...
                testNum,nTests,k(i),r(i));
    result(testNum) = poissonMMS_test(dom,params,options);

    % display result
    if result(testNum) == 1
        fprintf('Result: PASS\n\n');
    else
        fprintf('Result: FAIL\n\n');
    end

    % prepare for next test
    testNum = testNum + 1;

end


% SUMMARIZE RESULTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
executionTime = toc(start);

% gather all strings needed for consistent spacing and find max length
allStrs = cell(loopsize(1),1);
for i = 1:loopsize(1)
    allStrs{i} = sprintf(' coefficients: k = %s, r = %s.', k(i), r(i));
end
maxLen = max(cellfun(@length, allStrs));

% Print summary 
fprintf('SUMMARY for poissonMMS boundary condition test\n');
for i = 1:loopsize(1)
    if result(i) == 1
        passFail = 'PASS';
    else
        passFail = 'FAIL';
    end
    
    % %-*s formats the string to a fixed width based on maxLen
    fprintf('%-*s   Result: %s\n', maxLen, allStrs{i}, passFail);
end

fprintf('Total execution time: %f s\n\n', executionTime);
