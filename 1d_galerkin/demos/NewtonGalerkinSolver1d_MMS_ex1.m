% NOTE: DO NOT CHANGE THE SOLUTION OR PARAMETERS
% If you wish to test other problems, create a new script based on this one.
%
% This script demonstrates the failure of the NewtonGalerkinSolver1d_poisson
% class to converge for a particular manufactured solution. 
% 
% SETUP:
% The problem solved is:
%
%    -d/dx(k du/dx) + r u = f,  x in (0,1)
%                -k du/dn = g_1, x = 0
%                -k du/dn = g_2, x = 1
%
% Here, g_1 and g_2 can either be manufactured Neumann boundary conditions, i.e.
% linear boundary conditions where the values g_1 and g_2 are determined from
% the desired solution uTrue. Or g_1 and g_2 can be nonlinear boundary condtions
% given by:
% 
%    g_1(u) = exp(u)
%    g_2(u) = -(exp(u-2) + 2)
%
% These nonlinear boundary conditions have been chosen so to be compatible with
% the true solution
% 
%    uTrue = x + x^2.
%
% The initial guess is set to be either the constant function u = 1, or a
% perturbation of the true solution.
%
%
% RESULTS:
% The choice of:
%
%    LEFT_nonlinearBCs = 0
%    RIGHT_nonlinearBCs = 1
%    constant_initial_guess = 1
%    rU_linear = 1
% 
% will cause the Newton iteration to fail to converge. The choice of 
%
%   LEFT_nonlinearBCs = 1 
%   RIGHT_nonlinearBCs = 1
%   constant_initial_guess = 1
%   rU_linear = 1
%
% will cause order 1 convergence. All other choices cause order 2 convergence.
%
% EXPLANATION:
% The issue appears to be that the initial guess u = 1 is outside the basin of
% convergence for Newton's method when the BC on the right side is nonlinear,
% the BC on the left side is linear, and r(u) = u. Changing any of these three
% conditions appears to place the initial guess within the basin of convergence.
%
% Intiutively, this makes some sense. When r(u) = u, the problem is linear on
% the interior and the nonlinearity is on only the right boundary node. The
% nonlinear boundary condition there is exponential. So if the initial guess is
% too far from the true solution, then Newton's method actually converges to a
% different (and spurious) solution. If both boundary conditions are nonlinear,
% then the nonlinearity on the left boundary appears to help guide the solution
% towards the true solution. If r(u) = u^2, then the problem is nonlinear on the
% interior as well, which also appears to help guide the solution towards the
% true solution.


clear all; syms x t u;
% USER INPUTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% original problem
uTrue = x + x^2;

% determine to use nonlinear BCs on left, right, or both
LEFT_nonlinearBCs = 1;      % if 1, use nonlinear BCs; if 0, use linear (manufactured) BCs
RIGHT_nonlinearBCs = 1;     % if 1, use nonlinear BCs; if 0, use linear (manufactured) BCs
constant_initial_guess = 1; % if 1, use constant initial guess u = 1; else use perturbed true solution
rU_linear = 1;              % if 1, use r(u) = u; if 0, use r(u) = u^2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% DOMAIN PARAMETERS
% domain parameters
xLim = [0 1];

% mesh parameters
base = 2;

% PDE PARAMETERS
% pde coefficients
k = 1;
if rU_linear
    r = u;
else
    r = u^2;
end

% boundary conditions
BCtypes = 'NN';

% original problem
U_N_L = exp(u);
dU_N_L = exp(u);
U_N_R = -(exp(u-2) + 2);
dU_N_R = -exp(u-2);

BCvals = {U_N_L, U_N_R};
BCvals_du = {dU_N_L, dU_N_R};

% BLACK BOX %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% create domain 
dom = Domain1d(xLim);

% manufacture auxiliary functions
manfun = ManufacturedFunctions1d_poisson(k,r,uTrue,BCtypes);
auxfun = manfun.outputCoefficients();

% manufacture boundary conditions, then replace one or both with nonlinear conditions
%   Note: manfun.outputBoundaryConditions creates linear (manufactured) BCs. These
%   can then be manually overwritten on the left or right edges. The setter for
%   the BoundaryEdge1d class can convert symbolic expressions into function
%   handles as needed.
dom.boundary = manfun.outputBoundaryConditions(dom,BCtypes);

if LEFT_nonlinearBCs
    dom.boundary.edges(1).boundaryCondition = U_N_L;
    dom.boundary.edges(1).boundaryCondition_ddu = dU_N_L;
end

if RIGHT_nonlinearBCs
    dom.boundary.edges(2).boundaryCondition = U_N_R;
    dom.boundary.edges(2).boundaryCondition_ddu = dU_N_R;
end

%{
% create small problem for checking matrices
dom.mesh = Mesh1d(xLim,2,1);
prob = NewtonGalerkinSolver1d_poisson(dom,auxfun);
uTrue_handle = matlabFunction(symfun(uTrue,x));
prob.U = uTrue_handle(prob.domain.mesh.nodes); % change initial guess to be the exact solution
prob = prob.assembleTensors;
prob = prob.assembleVectors;
prob = prob.assembleBCs;
[DJ,J] = prob.finalAssembly(prob.U);
%prob = prob.solve;

return
%}

fprintf('Begin MMS test\n');
% loop over mesh refinements
err = {};
i = 0;
for p = 2:5

    fprintf('  p = %d:',p); tic;
    % update counter
    i = i + 1;

    % update mesh
    dom.mesh = Mesh1d(xLim,base,p);

    % create galerkin solver
    prob = NewtonGalerkinSolver1d_poisson(dom,auxfun);

    % manually set initial guess
    %   NOTE: the initial guess is automatically set to U = 1 in the
    %   NewtonGalerkinSolver1d class. But it can be overwritten using the
    %   following code.
    if ~constant_initial_guess
        uTrue_handle = matlabFunction(symfun(uTrue,x));
        prob.U = uTrue_handle(prob.domain.mesh.nodes);
        prob.U = prob.U + 0.1; % perturb initial guess slightly
    end

    % solve
    prob = prob.solve;

    % compute and display errors
    err{i} = dom.L2err_twoPointQuadrature(prob.solution,uTrue);

    % compute difference in errors
    if i > 1
        rate{i-1} = log(err{i-1}/err{i})/log(base);
    end
    timeRequired = toc;
    fprintf(' %.4f s\n',timeRequired);

    prob.plot;
    pause()

end

% display results as table using matlab's table function
pList = (2:5).';
errList = cell2mat(err.');
rateList = [NaN; cell2mat(rate.')];
T = table(pList, errList, rateList);
T.Properties.VariableNames = {'p' 'L2_Error' 'Convergence_Rate'};
disp(T);
