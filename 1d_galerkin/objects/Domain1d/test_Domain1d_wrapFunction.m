% test_Domain1d_wrapFunction.m
%   Test script for Domain1d.function2nodes method
%   - tests constant, space-varying, time-varying, and space-and-time-varying functions
%   - tests different input function types: double, function_handle, symfun, sym
%   - reports pass/fail for each test case (1 = pass, 0 = fail; but
%           realistically it will just throw an error on fail)
%
% NOTES:
%   (1) Function evaluations in Domain1d assume that f is a function handle.
%   However, this package allows you to pass functions as doubles, syms, or
%   symfuns as well. In these cases, the function is converted to a function
%   handle internally before evaluation. Therefore, this test script includes
%   tests for these other input types as well.
%   
%   (2) Function evaluations assume that f = f(x,t). However, if f = f(), f =
%   f(x), or f = f(t), these cases are handled internally by converting to f =
%   f(x,t) before evaluation. Therefore, this test script includes tests for
%   these other input variable cases as well.
%   
%   (3) Some combinations of input function types and variable dependencies are
%   not possible to implement. In these cases, the test is automatically marked
%   as passed. For example:
%       - A double can only represent a constant function. 
%       - A symfun must have some variable defined. So it is not possible for a
%           symfun to represent a constant function with no variables.  
%       - A sym is defined using like, e.g. f = sin(2*pi*x). So it is not
%           possible for this sym, which is a function of space, to have
%           variables of space and time.  
%   We still note these cases in the output for completeness.
%
%   (4) The expected output for a function evaluation in this package is an
%   array of size nSpaceEvals x nTimeEvals. So we check that the output size is
%   correct as well as the output values. This is especially important for
%   functions that are constant in either space or time, since MATLAB's default
%   behavior is to return a single value in these cases rather than an array of
%   repeated values. This behavior is addressed in the Domain1d.wrapFunctions
%   method by adding an array of zeros of the correct size to the output.
%
%   (5) This test script also includes a test for a Domain1d with no time
%   stepping (i.e., static problem). Internally, this case is handled by
%   assuming that t = 0, meaning that nTimeEvals = 1, and therefore the result
%   will be a vector of size nSpaceEvals x 1. Therefore, if a time-varying
%   function is passed to a Domain1d with no time stepping, the function is
%   evaluated at t = 0.

clear all; syms x t;
test_constants          = 1;
test_spaceVarying       = 1;
test_timeVarying        = 1;
test_spaceTimeVarying   = 1;
test_noTimeStepping     = 1;

% USER INPUTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% domain parameters
xLim = [0 1];

% mesh parameters
base = 2;
p = 5;

% time stepping parameters
T = 1;
dt = 0.1;


% BLACK BOX %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% create domain and mesh
dom = Domain1d(xLim);
mesh = Mesh1d(xLim,base,p);
time = TimeStepping(T,dt);
dom.mesh = mesh;
dom.time = time;

if test_constants

    fprintf('--- Testing constant functions ---\n');

    % TEST CONSTANT FUNCTIONS
    test_val = 10 * ones(dom.mesh.nNodes,dom.time.N_t);

    %   double
    func = 10;
    pass = min(min(dom.function2nodes(func) == test_val));
    fprintf('%d: Constant function, double\n', pass);

    %   function_handle, no vars
    func = @( ) 10;
    pass = min(min(dom.function2nodes(func) == test_val));
    fprintf('%d: Constant function, function_handle no vars\n', pass);

    %   function_handle, x var
    func = @(x) 10;
    pass = min(min(dom.function2nodes(func) == test_val));
    fprintf('%d: Constant function, function_handle x var\n', pass);

    %   function_handle, t var
    func = @(t) 10;
    pass = min(min(dom.function2nodes(func) == test_val));   
    fprintf('%d: Constant function, function_handle t var\n', pass);

    %   function_handle, x and t vars
    func = @(x,t) 10;
    pass = min(min(dom.function2nodes(func) == test_val));
    fprintf('%d: Constant function, function_handle x and t vars\n', pass);

    %   symfun, no vars
    fprintf('%d: Constant function, symfun no vars: NOT POSSIBLE (automatic pass)\n',1);

    %   symfun, x var
    func = symfun(10, x);
    pass_size = min(min(size(dom.function2nodes(func)) == size(test_val)));
    pass_val = min(min(dom.function2nodes(func) == test_val));
    fprintf('%d: Constant function, symfun x vard\n', min(pass_size, pass_val));

    %   symfun, t var
    func = symfun(10, t);
    pass_size = min(min(size(dom.function2nodes(func)) == size(test_val)));
    pass_val = min(min(dom.function2nodes(func) == test_val));
    fprintf('%d: Constant function, symfun t var\n', min(pass_size, pass_val));

    %   symfun, x and t vars
    func = symfun(10, [x, t]);
    pass_size = min(min(size(dom.function2nodes(func)) == size(test_val)));
    pass_val = min(min(dom.function2nodes(func) == test_val));
    fprintf('%d: Constant function, symfun x and t vars\n', min(pass_size, pass_val));

    %   sym, no vars
    func = sym(10);
    pass_size = min(min(size(dom.function2nodes(func)) == size(test_val)));
    pass_val = min(min(dom.function2nodes(func) == test_val));
    fprintf('%d: Constant function, sym no vars\n', min(pass_size, pass_val));

    %   sym, x var
    fprintf('%d: Constant function, sym x var: NOT POSSIBLE (automatic pass)\n',1);

    %   sym, t var
    fprintf('%d: Constant function, sym t var: NOT POSSIBLE (automatic pass)\n',1);

    %   sym, x and t vars
    fprintf('%d: Constant function, sym x and t vars: NOT POSSIBLE (automatic pass)\n',1);
    fprintf('\n');
end


if test_spaceVarying

    fprintf('--- Testing space-varying functions ---\n');

    % TEST SPATIAL FUNCTIONS
    test_val = repmat(sin(2*pi*dom.mesh.nodes),1,dom.time.N_t);

    %   double
    fprintf('%d: Function of space, double: NOT POSSIBLE (automatic pass)\n',1);

    %   function_handle, no vars
    fprintf('%d: Function of space, function_handle no vars: NOT POSSIBLE (automatic pass)\n',1);

    %   function_handle, x var
    func = @(x) sin(2*pi*x);
    pass_size = min(min(size(dom.function2nodes(func)) == size(test_val)));
    pass_val = min(min(dom.function2nodes(func) == test_val));
    fprintf('%d: Function of space, function_handle x var: %d\n', min(pass_size, pass_val));

    %   function_handle, t var
    fprintf('%d: Function of space, function_handle t var: NOT POSSIBLE (automatic pass)\n',1);

    %   function_handle, x and t vars
    func = @(x,t) sin(2*pi*x);
    pass_size = min(min(size(dom.function2nodes(func)) == size(test_val)));
    pass_val = min(min(dom.function2nodes(func) == test_val));
    fprintf('%d: Function of space, function_handle x and t vars\n', min(pass_size, pass_val)); 

    %   symfun, no vars
    fprintf('%d: Function of space, symfun no vars: NOT POSSIBLE (automatic pass)\n',1);

    %   symfun, x var
    func = symfun(sin(2*pi*x), x);
    pass_size = min(min(size(dom.function2nodes(func)) == size(test_val)));
    pass_val = min(min(dom.function2nodes(func) == test_val));
    fprintf('%d: Function of space, symfun x var\n', min(pass_size, pass_val));

    %   symfun, t var
    fprintf('%d: Function of space, symfun t var: NOT POSSIBLE (automatic pass)\n',1);
    %   symfun, x and t vars
    func = symfun(sin(2*pi*x), [x, t]);
    pass_size = min(min(size(dom.function2nodes(func)) == size(test_val)));
    pass_val = min(min(dom.function2nodes(func) == test_val));
    fprintf('%d: Function of space, symfun x and t vars\n', min(pass_size, pass_val));

    %   sym, no vars
    fprintf('%d: Function of space, sym no vars: NOT POSSIBLE (automatic pass)\n',1);

    %   sym, x var
    func = sin(2*pi*x);
    pass_size = min(min(size(dom.function2nodes(func)) == size(test_val)));
    pass_val = min(min(dom.function2nodes(func) == test_val));
    fprintf('%d: Function of space, sym x var\n', min(pass_size, pass_val));

    %   sym, t var
    fprintf('%d: Function of space, sym t var: NOT POSSIBLE (automatic pass)\n',1);

    %   sym, x and t vars
    fprintf('%d: Function of space, sym x and t vars: NOT POSSIBLE (automatic pass)\n',1);
    fprintf('\n');
end


if test_timeVarying

    fprintf('--- Testing time-varying functions ---\n');

    % TEST FUNCTIONS OF TIME ONLY
    test_val = repmat(sin(2*pi*dom.time.tGrid),dom.mesh.nNodes,1);

    %   double
    fprintf('%d: Function of time, double: NOT POSSIBLE (automatic pass)\n',1);

    %   function_handle, no vars
    fprintf('%d: Function of time, function_handle no vars: NOT POSSIBLE (automatic pass)\n',1);

    %   function_handle, x var
    fprintf('%d: Function of time, function_handle x var: NOT POSSIBLE (automatic pass)\n',1);

    %   function_handle, t var
    func = @(t) sin(2*pi*t);
    pass_size = min(min(size(dom.function2nodes(func)) == size(test_val)));
    pass_val = min(min(dom.function2nodes(func) == test_val));
    fprintf('%d: Function of time, function_handle t var\n', min(pass_size, pass_val));

    %   function_handle, x and t vars
    func = @(x,t) sin(2*pi*t);
    pass_size = min(min(size(dom.function2nodes(func)) == size(test_val)));
    pass_val = min(min(dom.function2nodes(func) == test_val));
    fprintf('%d: Function of time, function_handle x and t vars\n', min(pass_size, pass_val));

    %   symfun, no vars
    fprintf('%d: Function of time, symfun no vars: NOT POSSIBLE (automatic pass)\n',1);

    %   symfun, x var
    fprintf('%d: Function of time, symfun x var: NOT POSSIBLE (automatic pass)\n',1);

    %   symfun, t var
    func = symfun(sin(2*pi*t), t);
    pass_size = min(min(size(dom.function2nodes(func)) == size(test_val)));
    pass_val = min(min(dom.function2nodes(func) == test_val));
    fprintf('%d: Function of time, symfun t var\n', min(pass_size, pass_val));

    %   symfun, x and t vars
    func = symfun(sin(2*pi*t), [x, t]);
    pass_size = min(min(size(dom.function2nodes(func)) == size(test_val)));
    pass_val = min(min(dom.function2nodes(func) == test_val));
    fprintf('%d: Function of time, symfun x and t vars\n', min(pass_size, pass_val));   

    %   sym, x var
    fprintf('%d: Function of time, sym x var: NOT POSSIBLE (automatic pass)\n',1);

    %   sym, t var
    func = sin(2*pi*t);
    pass_size = min(min(size(dom.function2nodes(func)) == size(test_val)));
    pass_val = min(min(dom.function2nodes(func) == test_val));
    fprintf('%d: Function of time, sym t var\n', min(pass_size, pass_val));

    %   sym, x and t vars
    fprintf('%d: Function of time, sym x and t vars: NOT POSSIBLE (automatic pass)\n',1);

    fprintf('\n');
end



if test_spaceTimeVarying

    fprintf('--- Testing space-and-time-varying functions ---\n');

    % TEST FUNCTIONS OF TIME AND SPACE
    test_val = sin(2*pi*dom.mesh.nodes) * exp(-dom.time.tGrid);

    %   double
    fprintf('%d: Function of space and time, double: NOT POSSIBLE (automatic pass)\n',1);

    %   function_handle, no vars
    fprintf('%d: Function of space and time, function_handle no vars: NOT POSSIBLE (automatic pass)\n',1);

    %   function_handle, x var
    fprintf('%d: Function of space and time, function_handle x var: NOT POSSIBLE (automatic pass)\n',1);

    %   function_handle, t var
    fprintf('%d: Function of space and time, function_handle t var: NOT POSSIBLE (automatic pass)\n',1);

    %   function_handle, x and t vars
    func = @(x,t) sin(2*pi*x) * exp(-t);
    pass_size = min(min(size(dom.function2nodes(func)) == size(test_val)));
    pass_val = min(min(dom.function2nodes(func) == test_val));
    fprintf('%d: Function of space and time, function_handle x and t vars\n', min(pass_size, pass_val));

    %   symfun, no vars
    fprintf('%d: Function of space and time, symfun no vars: NOT POSSIBLE (automatic pass)\n',1);

    %   symfun, x var
    fprintf('%d: Function of space and time, symfun x var: NOT POSSIBLE (automatic pass)\n',1);

    %   symfun, t var
    fprintf('%d: Function of space and time, symfun t var: NOT POSSIBLE (automatic pass)\n',1);

    %   symfun, x and t vars
    func = symfun(sin(2*pi*x) * exp(-t), [x, t]);
    pass_size = min(min(size(dom.function2nodes(func)) == size(test_val)));
    pass_val = min(min(dom.function2nodes(func) == test_val));
    fprintf('%d: Function of space and time, symfun x and t vars\n', min(pass_size, pass_val));

    %   sym, no vars
    fprintf('%d: Function of space and time, sym no vars: NOT POSSIBLE (automatic pass)\n',1);

    %   sym, x var
    fprintf('%d: Function of space and time, sym x var: NOT POSSIBLE (automatic pass)\n',1);

    %   sym, t var
    fprintf('%d: Function of space and time, sym t var: NOT POSSIBLE (automatic pass)\n',1);

    %   sym, x and t vars
    func = sin(2*pi*x) * exp(-t);
    pass_size = min(min(size(dom.function2nodes(func)) == size(test_val)));
    pass_val = min(min(dom.function2nodes(func) == test_val));
    fprintf('%d: Function of space and time, sym x and t vars\n', min(pass_size, pass_val));
    fprintf('\n');
end

if test_noTimeStepping

    fprintf('--- Testing domain with no time stepping ---\n');

    % create domain and mesh
    dom = Domain1d(xLim);
    mesh = Mesh1d(xLim,base,p);
    dom.mesh = mesh;

    % TEST FUNCTIONS OF TIME AND SPACE
    test_val = sin(2*pi*dom.mesh.nodes);

    %   double
    test_val = 10 * ones(dom.mesh.nNodes,1);
    func = 10;
    pass = min(dom.function2nodes(func) == test_val);
    fprintf('%d: Constant function, double\n', pass);

    %   function_handle, no vars
    func = @( ) 10;
    pass = min(dom.function2nodes(func) == test_val);
    fprintf('%d: Constant function, function_handle no vars\n', pass);

    %   function_handle, x var
    test_val = sin(2*pi*dom.mesh.nodes);
    func = @(x) sin(2*pi*x);
    pass = min(dom.function2nodes(func) == test_val);
    fprintf('%d: Function of space, function_handle x var\n', pass);
    
    % function_handl, t var
    test_val = 10 * exp(-0) * ones(dom.mesh.nNodes,1);
    func = @(t) 10 * exp(-t);
    pass = min(dom.function2nodes(func) == test_val);
    fprintf('%d: Function of time, function_handle t var\n', pass);

    %   function_handle, x and t vars
    test_val = sin(2*pi*dom.mesh.nodes) * exp(-0);
    func = @(x,t) sin(2*pi*x) * exp(-t);
    pass = min(dom.function2nodes(func) == test_val);
    fprintf('%d: Function of space and time, function_handle x and t vars\n', pass);

    %   symfun, no vars
    fprintf('%d: Constant function, symfun no vars: NOT POSSIBLE (automatic pass)\n',1);

    %   symfun, x var
    test_val = sin(2*pi*dom.mesh.nodes);
    func = symfun(sin(2*pi*x), x);
    pass = min(dom.function2nodes(func) == test_val);
    fprintf('%d: Function of space, symfun x var\n', pass);

    %   symfun, t var
    test_val = 10 * exp(-0) * ones(dom.mesh.nNodes,1);
    func = symfun(10 * exp(-t), t);
    pass = min(dom.function2nodes(func) == test_val);
    fprintf('%d: Function of time, symfun t var\n', pass);

    %   symfun, x and t vars
    test_val = sin(2*pi*dom.mesh.nodes) * exp(-0);
    func = symfun(sin(2*pi*x) * exp(-t), [x, t]);
    pass = min(dom.function2nodes(func) == test_val);
    fprintf('%d: Function of space and time, symfun x and t vars\n', pass);

    %   sym, no vars
    fprintf('%d: Constant function, sym no vars: NOT POSSIBLE (automatic pass)\n',1);

    %   sym, x var
    test_val = sin(2*pi*dom.mesh.nodes);
    func = sin(2*pi*x);
    pass = min(dom.function2nodes(func) == test_val);
    fprintf('%d: Function of space, sym x var\n', pass);

    %   sym, t var
    test_val = 10 * exp(-0) * ones(dom.mesh.nNodes,1);
    func = 10 * exp(-t);
    pass = min(dom.function2nodes(func) == test_val);
    fprintf('%d: Function of time, sym t var\n', pass); 

    %   sym, x and t vars
    test_val = sin(2*pi*dom.mesh.nodes) * exp(-0);
    func = sin(2*pi*x) * exp(-t);
    pass = min(dom.function2nodes(func) == test_val);
    fprintf('%d: Function of space and time, sym x and t vars\n', pass);

    fprintf('\n');
end