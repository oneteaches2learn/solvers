% rxndiffMMS_test_BCs
% Cycles through different domain types
clear all; x = sym('x',[1 2],'real'); syms t; start = tic;

% which tests to run:
inclusionTypes   = 1;
inclusionNumbers = 1;
inclusionRatios  = 1;
oblongDomains    = 1;
yLines           = 1;
effectiveRegions = 1;

% FIXED PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% domain bounds
% params.xLim_dom = ...
% params.yLim_dom = ...
% ... set in loop for this test

% Y bounds
% params.xLim_Y = ...
% params.yLim_Y = ...
% ... set in loop for this test

% number of inclusions
% params.incRatio = ...
% params.eps = ...
% ... set in loop for this test

% time stepping
params.T = 1;

% boundary conditions
params.bTypes_outer = 'RRRR';
params.bTypes_inner = 'R';

% coefficients
params.c = 2 + sin(x(1) + x(2)) * t^2 + t^2;
params.k = 2 + sin(x(1) + x(2)) * t^2 + t^2;
params.r = 2 + sin(x(1) + x(2)) * t^2 + t^2;

% desired true solution
params.uTrue = cos(2 * pi * x(1)) * cos(2 * pi * x(2));

% mms parameters
optionSet.base = 2;
optionSet.demo = 0;
optionSet.timeOffset = 4;
optionSet.timeFactor = 2;
optionSet.pmin = 4;
optionSet.pmax = 7;
optionSet.meshInclusions = 'on';
optionSet.effectiveRegion = 'Omega_eps';
optionSet.expectedRate = 2;
optionSet.convergenceTolerance = 0.1;


% EXPIREMENTAL PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% coefficients
% create domain
fprintf('Initialization\n')
fprintf(' Contructing Domains:'), tic
idx = 1;

% TESTS OF INCLUSION TYPES
if inclusionTypes == 1
    params.xLim_dom = [0, 1];
    params.yLim_dom = [0, 1];
    params.xLim_Y = [0, 1];
    params.yLim_Y = [0, 1];
    params.incRatio = 1;
    params.eps = 1;

    % test 1: No inclusions
    str{idx} = 'No inclusions';
    dom = Domain2d(params.xLim_dom,params.yLim_dom);
    dom = dom.setBCTypes([params.bTypes_outer,params.bTypes_inner]);
    dom = dom.setTime(params.T);
    domains{idx} = dom;
    options{idx} = optionSet;
    idx = idx + 1;

    % test 2: Circle inclusion
    str{idx} = 'One circular inclusion';
    inc = Inclusion2d_circle(params.xLim_Y,params.yLim_Y,params.incRatio);
    dom = Domain2d_punctured(params.xLim_dom,params.yLim_dom,inc,params.eps);
    dom = dom.setBCTypes([params.bTypes_outer,params.bTypes_inner]);
    dom = dom.setTime(params.T);
    domains{idx} = dom;
    options{idx} = optionSet;
    idx = idx + 1;

    % test 3: Square inclusion
    str{idx} = 'One square inclusion';
    inc = Inclusion2d_square(params.xLim_Y,params.yLim_Y,params.incRatio);
    dom = Domain2d_punctured(params.xLim_dom,params.yLim_dom,inc,params.eps);
    dom = dom.setBCTypes([params.bTypes_outer,params.bTypes_inner]);
    dom = dom.setTime(params.T);
    domains{idx} = dom;
    options{idx} = optionSet;
    idx = idx + 1;
end


% TESTS OF INCLUSION NUMBERS
if inclusionNumbers == 1
    params.xLim_dom = [0, 1];
    params.yLim_dom = [0, 1];
    params.xLim_Y = [0, 1];
    params.yLim_Y = [0, 1];
    params.incRatio = 1;
    params.eps = 1/3;

    % test 4: Nine circular inclusions
    str{idx} = 'Nine circular inclusions';
    inc = Inclusion2d_circle(params.xLim_Y,params.yLim_Y,params.incRatio);
    dom = Domain2d_punctured(params.xLim_dom,params.yLim_dom,inc,params.eps);
    dom = dom.setBCTypes([params.bTypes_outer,params.bTypes_inner]);
    dom = dom.setTime(params.T);
    domains{idx} = dom;
    options{idx} = optionSet;
    idx = idx + 1;

    % test 5: Nine Square inclusions
    str{idx} = 'Nine square inclusions';
    inc = Inclusion2d_square(params.xLim_Y,params.yLim_Y,params.incRatio);
    dom = Domain2d_punctured(params.xLim_dom,params.yLim_dom,inc,params.eps);
    dom = dom.setBCTypes([params.bTypes_outer,params.bTypes_inner]);
    dom = dom.setTime(params.T);
    domains{idx} = dom;
    options{idx} = optionSet;
    idx = idx + 1;
end


% TESTS OF INCLUSION RATIOS
if inclusionRatios == 1
    params.xLim_dom = [0, 2];
    params.yLim_dom = [0, 1];
    params.xLim_Y = [0, 1];
    params.yLim_Y = [0, 1];
    params.incRatio = 1;
    params.eps = 1;

    % test 6: Small circular inclusion
    str{idx} = 'Small circular inclusion (incRatio = 1/3)';
    params.incRatio = 1/3;
    inc = Inclusion2d_circle(params.xLim_Y,params.yLim_Y,params.incRatio);
    dom = Domain2d_punctured(params.xLim_dom,params.yLim_dom,inc,params.eps);
    dom = dom.setBCTypes([params.bTypes_outer,params.bTypes_inner]);
    dom = dom.setTime(params.T);
    domains{idx} = dom;
    options{idx} = optionSet;
    idx = idx + 1;

    % test 7: Small square inclusion
    str{idx} = 'Small circular inclusion (incRatio = 1/3)';
    params.incRatio = 1/3;
    inc = Inclusion2d_square(params.xLim_Y,params.yLim_Y,params.incRatio);
    dom = Domain2d_punctured(params.xLim_dom,params.yLim_dom,inc,params.eps);
    dom = dom.setBCTypes([params.bTypes_outer,params.bTypes_inner]);
    dom = dom.setTime(params.T);
    domains{idx} = dom;
    options{idx} = optionSet;
    idx = idx + 1;

    % test 8: Large circular inclusion
    str{idx} = 'Large circular inclusion (incRatio = 3/2)';
    params.incRatio = 3/2;
    inc = Inclusion2d_circle(params.xLim_Y,params.yLim_Y,params.incRatio);
    dom = Domain2d_punctured(params.xLim_dom,params.yLim_dom,inc,params.eps);
    dom = dom.setBCTypes([params.bTypes_outer,params.bTypes_inner]);
    dom = dom.setTime(params.T);
    domains{idx} = dom;
    options{idx} = optionSet;
    idx = idx + 1;

    % test 9: Large square inclusion
    str{idx} = 'Large circular inclusion (incRatio = 3/2)';
    params.incRatio = 3/2;
    inc = Inclusion2d_square(params.xLim_Y,params.yLim_Y,params.incRatio);
    dom = Domain2d_punctured(params.xLim_dom,params.yLim_dom,inc,params.eps);
    dom = dom.setBCTypes([params.bTypes_outer,params.bTypes_inner]);
    dom = dom.setTime(params.T);
    domains{idx} = dom;
    options{idx} = optionSet;
    idx = idx + 1;
end


% TESTS OF OBLONG DOMAIN / INCLUSION SIZES
if oblongDomains == 1
    params.incRatio = 1;
    params.eps = 1;

    % test 10: Rectangular domain 1
    str{idx} = 'Rectangular domain (x wider dimension)';
    params.xLim_dom = [0, 2];
    params.yLim_dom = [0, 1];
    params.xLim_Y = [0, 1];
    params.yLim_Y = [0, 1];
    inc = Inclusion2d_circle(params.xLim_Y,params.yLim_Y,params.incRatio);
    dom = Domain2d_punctured(params.xLim_dom,params.yLim_dom,inc,params.eps);
    dom = dom.setBCTypes([params.bTypes_outer,params.bTypes_inner]);
    dom = dom.setTime(params.T);
    domains{idx} = dom;
    options{idx} = optionSet;
    idx = idx + 1;

    % test 11: Rectangular domain 2
    str{idx} = 'Rectangular domain (y wider dimension)';
    params.xLim_dom = [0, 1];
    params.yLim_dom = [0, 2];
    params.xLim_Y = [0, 1];
    params.yLim_Y = [0, 1];
    inc = Inclusion2d_circle(params.xLim_Y,params.yLim_Y,params.incRatio);
    dom = Domain2d_punctured(params.xLim_dom,params.yLim_dom,inc,params.eps);
    dom = dom.setBCTypes([params.bTypes_outer,params.bTypes_inner]);
    dom = dom.setTime(params.T);
    domains{idx} = dom;
    options{idx} = optionSet;
    idx = idx + 1;

    % test 12: Rectangular domain 3
    str{idx} = 'Rectangular Y (x wider dimension)';
    params.xLim_dom = [0, 1];
    params.yLim_dom = [0, 1];
    params.xLim_Y = [0, 1];
    params.yLim_Y = [0, 0.5];
    inc = Inclusion2d_circle(params.xLim_Y,params.yLim_Y,params.incRatio);
    dom = Domain2d_punctured(params.xLim_dom,params.yLim_dom,inc,params.eps);
    dom = dom.setBCTypes([params.bTypes_outer,params.bTypes_inner]);
    dom = dom.setTime(params.T);
    domains{idx} = dom;
    options{idx} = optionSet;
    idx = idx + 1;

    % test 13: Rectangular domain 4
    str{idx} = 'Rectangular Y (y wider dimension)';
    params.xLim_dom = [0, 1];
    params.yLim_dom = [0, 1];
    params.xLim_Y = [0, 0.5];
    params.yLim_Y = [0, 1];
    inc = Inclusion2d_circle(params.xLim_Y,params.yLim_Y,params.incRatio);
    dom = Domain2d_punctured(params.xLim_dom,params.yLim_dom,inc,params.eps);
    dom = dom.setBCTypes([params.bTypes_outer,params.bTypes_inner]);
    dom = dom.setTime(params.T);
    domains{idx} = dom;
    options{idx} = optionSet;
    idx = idx + 1;
end


% TESTS OF Y-LINE
if yLines == 1
    params.xLim_dom = [0, 1];
    params.yLim_dom = [0, 1];
    params.xLim_Y = [0, 1];
    params.yLim_Y = [0, 1];
    params.incRatio = 1;
    params.eps = 1/2;

    % test 14: One y-line 1
    str{idx} = 'One yLine (avoids inclusions)';
    inc = Inclusion2d_circle(params.xLim_Y,params.yLim_Y,params.incRatio);
    dom = Domain2d_punctured(params.xLim_dom,params.yLim_dom,inc,params.eps);
    dom = dom.add_yline(0.5);
    dom = dom.setBCTypes([params.bTypes_outer,params.bTypes_inner]);
    dom = dom.setTime(params.T);
    domains{idx} = dom;
    options{idx} = optionSet;
    idx = idx + 1;

    % test 15: One y-line 2
    str{idx} = 'One yLine (no inclusions)';
    dom = Domain2d(params.xLim_dom,params.yLim_dom);
    dom = dom.add_yline(0.5);
    dom = dom.setBCTypes([params.bTypes_outer,params.bTypes_inner]);
    dom = dom.setTime(params.T);
    domains{idx} = dom;
    options{idx} = optionSet;
    idx = idx + 1;

    % test 16: One y-line 3
    str{idx} = 'One yLine (intersects circular inclusion)';
    inc = Inclusion2d_circle(params.xLim_Y,params.yLim_Y,params.incRatio);
    dom = Domain2d_punctured(params.xLim_dom,params.yLim_dom,inc,params.eps);
    dom = dom.add_yline(0.75);
    dom = dom.setBCTypes([params.bTypes_outer,params.bTypes_inner]);
    dom = dom.setTime(params.T);
    domains{idx} = dom;
    options{idx} = optionSet;
    idx = idx + 1;

    % test 17: One y-line 4
    str{idx} = 'One yLine (intersects square inclusion)';
    inc = Inclusion2d_square(params.xLim_Y,params.yLim_Y,params.incRatio);
    dom = Domain2d_punctured(params.xLim_dom,params.yLim_dom,inc,params.eps);
    dom = dom.add_yline(0.75);
    dom = dom.setBCTypes([params.bTypes_outer,params.bTypes_inner]);
    dom = dom.setTime(params.T);
    domains{idx} = dom;
    options{idx} = optionSet;
    idx = idx + 1;

    % test 18: One y-line 5
    str{idx} = 'One yLine (tangent to circular inclusion)';
    params.incRatio = pi/2;
    params.eps = 1;
    inc = Inclusion2d_circle(params.xLim_Y,params.yLim_Y,params.incRatio);
    dom = Domain2d_punctured(params.xLim_dom,params.yLim_dom,inc,params.eps);
    dom = dom.add_yline(0.75);
    dom = dom.setBCTypes([params.bTypes_outer,params.bTypes_inner]);
    dom = dom.setTime(params.T);
    domains{idx} = dom;
    options{idx} = optionSet;
    idx = idx + 1;

    % test 19: One y-line 6
    str{idx} = 'One yLine (tangent to square inclusion)';
    params.incRatio = 2;
    params.eps = 1;
    inc = Inclusion2d_square(params.xLim_Y,params.yLim_Y,params.incRatio);
    dom = Domain2d_punctured(params.xLim_dom,params.yLim_dom,inc,params.eps);
    dom = dom.add_yline(0.75);
    dom = dom.setBCTypes([params.bTypes_outer,params.bTypes_inner]);
    dom = dom.setTime(params.T);
    domains{idx} = dom;
    options{idx} = optionSet;
    idx = idx + 1;
end


% TESTS OF EFFECTIVE REGIONS
if effectiveRegions == 1
    params.xLim_dom = [0, 1];
    params.yLim_dom = [0, 1];
    params.xLim_Y = [0, 1];
    params.yLim_Y = [0, 1];
    params.incRatio = 2;
    params.eps = 1;

    % test 20: Effective regions 1
    str{idx} = 'Inclusion Mesh: On.  Effective Region: Omega_eps';
    optionSet.meshInclusions = 'on';
    optionSet.effectiveRegion = 'Omega_eps';
    inc = Inclusion2d_square(params.xLim_Y,params.yLim_Y,params.incRatio);
    dom = Domain2d_punctured(params.xLim_dom,params.yLim_dom,inc,params.eps);
    dom = dom.add_yline(0.6);
    dom = dom.setBCTypes([params.bTypes_outer,params.bTypes_inner]);
    dom = dom.setTime(params.T);
    domains{idx} = dom;
    options{idx} = optionSet;
    idx = idx + 1;

    % test 21: Effective regions 2
    str{idx} = 'Inclusion Mesh: Off. Effective Region: Omega_eps';
    optionSet.meshInclusions = 'off';
    optionSet.effectiveRegion = 'Omega_eps';
    inc = Inclusion2d_square(params.xLim_Y,params.yLim_Y,params.incRatio);
    dom = Domain2d_punctured(params.xLim_dom,params.yLim_dom,inc,params.eps);
    dom = dom.add_yline(0.6);
    dom = dom.setBCTypes([params.bTypes_outer,params.bTypes_inner]);
    dom = dom.setTime(params.T);
    domains{idx} = dom;
    options{idx} = optionSet;
    idx = idx + 1;

    % test 22: Effective regions 3
    str{idx} = 'Inclusion Mesh: On.  Effective Region: Omega';
    optionSet.meshInclusions = 'on';
    optionSet.effectiveRegion = 'Omega';
    inc = Inclusion2d_square(params.xLim_Y,params.yLim_Y,params.incRatio);
    dom = Domain2d_punctured(params.xLim_dom,params.yLim_dom,inc,params.eps);
    dom = dom.add_yline(0.6);
    dom = dom.setBCTypes([params.bTypes_outer,params.bTypes_inner]);
    dom = dom.setTime(params.T);
    domains{idx} = dom;
    options{idx} = optionSet;
    idx = idx + 1;

    % test 23: Effective regions 4
    str{idx} = 'Inclusion Mesh: Off. Effective Region: Omega';
    optionSet.meshInclusions = 'off';
    optionSet.effectiveRegion = 'Omega';
    inc = Inclusion2d_square(params.xLim_Y,params.yLim_Y,params.incRatio);
    dom = Domain2d_punctured(params.xLim_dom,params.yLim_dom,inc,params.eps);
    dom = dom.add_yline(0.6);
    dom = dom.setBCTypes([params.bTypes_outer,params.bTypes_inner]);
    dom = dom.setTime(params.T);
    domains{idx} = dom;
    options{idx} = optionSet;
    idx = idx + 1;
end


executionTime = toc; 
fprintf(' %f s\n',executionTime)


% RUN TEST %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% generate test parameters
testNum = 1;
loopsize(1) = length(domains);
nTests = prod(loopsize);

for i = 1:loopsize(1), 

    % store experimental variable in params
    dom_i = domains{i};
    opt_i = options{i};

    % run test
    fprintf('TEST %d of %d. rxndiffMMS. domain: %s\n', ...
                testNum,nTests,str{i});
    result(testNum) = rxndiffMMS_test(dom_i,params,opt_i);

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
    allStrs{i} = sprintf(' domain: %s', str{i});
end
maxLen = max(cellfun(@length, allStrs));

% Print summary 
fprintf('SUMMARY for rxndiffMMS domain test\n');
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
