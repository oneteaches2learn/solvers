% pennesMMS_test_BCs
% Cycles through all boundary conditions
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
% params.bTypes_outer = ...
% params.bTypes_inner = ...
% ... set in loop for this test

% coefficients
params.c = 2 + sin(x(1) + x(2)) * t + t^2; 
params.k = 2 + sin(x(1) + x(2)) * t + t^2;
params.r = 2 + sin(x(1) + x(2)) * t + t^2;
params.uStar = 2 + sin(x(1) + x(2)) * t + t^2;

% desired true solution
params.uTrue = sin(2 * pi * x(1)) * sin(2 * pi * x(2)) * (t + 1) + (t + 1);

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
% outer boundary conditions
bTypes_outer(1,:) = ['DDDD'];
bTypes_outer(2,:) = ['NNNN'];
bTypes_outer(3,:) = ['RRRR'];

% inner boundary conditions
bTypes_inner(1,:) = ['D'];
bTypes_inner(2,:) = ['N'];
bTypes_inner(3,:) = ['R'];


% RUN TEST %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% generate test parameters
testNum = 1;
loopsize(1) = size(bTypes_outer,1);
loopsize(2) = size(bTypes_inner,1);
nTests = prod(loopsize);

for i = 1:loopsize(1), for j = 1:loopsize(2)

    % store experimental variable in params
    params.bTypes_outer = bTypes_outer(i,:);
    params.bTypes_inner = bTypes_inner(j,:);

    % create domain
    inc = Inclusion2d_circle(params.xLim_Y,params.yLim_Y,params.incRatio);
    %inc = Inclusion2d_square(params.xLim_Y,params.yLim_Y,params.incRatio);
    dom = Domain2d_punctured(params.xLim_dom,params.yLim_dom,inc,params.eps);
    %dom = Domain2d(xLim_dom,yLim_dom);
    dom = dom.add_yline(0.8);
    dom = dom.setBCTypes([params.bTypes_outer,params.bTypes_inner]);
    dom = dom.setTime(params.T);

    % run test
    fprintf('TEST %d of %d. pennesMMS. boundary conditions: %s, %s\n', ...
                testNum,nTests,bTypes_outer(i,:),bTypes_inner(j,:));
    result(testNum) = pennesMMS_test(dom,params,options);

    % display result
    if result(testNum) == 1
        fprintf('Result: PASS\n\n');
    else
        fprintf('Result: FAIL\n\n');
    end

    % prepare for next test
    testNum = testNum + 1;

end, end

% summarize results
executionTime = toc(start);
fprintf('SUMMARY for pennesMMS boundary condition test\n');
for i = 1:loopsize(1), for j = 1:loopsize(2)

    fprintf(' boundary conditions: %s, %s.\tResult: ',bTypes_outer(i,:),bTypes_inner(j,:));
    if result(i) == 1
        fprintf('PASS\n');
    else
        fprintf('FAIL\n');
    end
end, end
fprintf('Total execution time: %f s\n\n',executionTime);