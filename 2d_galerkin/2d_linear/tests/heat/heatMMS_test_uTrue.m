% rxndiffMMS_test_BCs
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
% params.bTypes_outer = ... set in loop for this test
params.bTypes_inner = 'R';

% coefficients
params.c = 2 + sin(x(1) + x(2)) + t^2;
params.k = 2 + sin(x(1) + x(2)) + t^2;

% desired true solution
% params.uTrue = ...
% ... set in loop for this test

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
% desired true solution
uTrue(1) = sin(2 * pi * x(1)) * sin(2 * pi * x(2));
uTrue(2) = cos(2 * pi * x(1)) * cos(2 * pi * x(2));
uTrue(3) = sin(2 * pi * x(1)) * sin(2 * pi * x(2)) * (t + 1) + (t + 1);
uTrue(4) = 0;
uTrue(5) = 1;
uTrue(6) = t;

% boundary conditions
bTypes_outer(1,:) = ['DDDD'];
bTypes_outer(2,:) = ['NNNN'];
bTypes_outer(3,:) = ['RRRR'];


% RUN TEST %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% generate test parameters
testNum = 1;
loopsize(1) = length(uTrue);
loopsize(2) = size(bTypes_outer,1);
nTests = prod(loopsize);

for i = 1:loopsize(1), for j = 1:loopsize(2)

    % store experimental variable in params
    params.uTrue = uTrue(i);
    params.bTypes_outer = bTypes_outer(j,:);

    % create domain
    inc = Inclusion2d_circle(params.xLim_Y,params.yLim_Y,params.incRatio);
    %inc = Inclusion2d_square(params.xLim_Y,params.yLim_Y,params.incRatio);
    dom = Domain2d_punctured(params.xLim_dom,params.yLim_dom,inc,params.eps);
    %dom = Domain2d(xLim_dom,yLim_dom);
    dom = dom.add_yline(0.8);
    dom = dom.setBCTypes([params.bTypes_outer,params.bTypes_inner]);
    dom = dom.setTime(params.T);

    % run test
    fprintf('TEST %d of %d. heatMMS. uTrue: %s, BCs: %s\n', ...
                testNum,nTests,uTrue(i),bTypes_outer(j,:));
    %{
    try
        result(testNum) = heatMMS_test(params,options);
    catch
        result(testNum) = 0;
    end
    %}
    result(testNum) = heatMMS_test(dom,params,options);

    % display result
    if result(testNum) == 1
        fprintf('Result: PASS\n\n');
    else
        fprintf('Result: FAIL\n\n');
    end

    % prepare for next test
    testNum = testNum + 1;

end, end


% SUMMARIZE RESULTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% gather all strings needed for consistent spacing in two columns
uTrueStr = cell(loopsize(1), loopsize(2));
BCsStr   = cell(loopsize(1), loopsize(2));

% Store each line's "uTrue:" text and "BCs:" text separately
for i = 1:loopsize(1)
    for j = 1:loopsize(2)
        uTrueStr{i, j} = sprintf('uTrue: %s', uTrue(i));
        BCsStr{i, j}   = sprintf('BCs: %s', bTypes_outer(j,:));
    end
end

% Find max width for each column
uTrueLen = max(cellfun(@length, uTrueStr(:)));
BCsLen   = max(cellfun(@length, BCsStr(:)));

% Print summary
fprintf('SUMMARY for heatMMS boundary condition test\n');

% We assume 'result' is a 1D array or can be indexed consistently here.
% If you store multiple results (one per (i,j)), adapt accordingly.
idx = 1;  % index into result
for i = 1:loopsize(1)
    for j = 1:loopsize(2)
        if result(idx) == 1
            passFail = 'PASS';
        else
            passFail = 'FAIL';
        end
        
        % %-*s aligns each column to the max width found above
        fprintf('%-*s   %-*s   Result: %s\n', ...
            uTrueLen, uTrueStr{i, j}, ...
            BCsLen,   BCsStr{i, j},   ...
            passFail);
        
        idx = idx + 1;
    end
end

executionTime = toc(start);
fprintf('Total execution time: %f s\n\n', executionTime);
