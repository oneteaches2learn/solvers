% For this test, we choose u_N_R = -exp(u-2) - 2 and uTrue = x + x^2. These
% choices "match" in the following way: 
% 
% -n cdot grad(uTrue) = - (1 + 2x) = -3 at x = 1, and
% - exp(uTrue(1)-2) - 2 = -exp(1+1-2) - 2 = -exp(0) - 2 = -1 - 2 = -3 at x = 1.
% 
% So the strategy is to choose a desired nonlinear boundary condition to test,
% then choose a uTrue that pointwise satisfies that boundary condition.

clear all; syms x t u;
% USER INPUTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% desired uTrue
uTrue = x + x^2;

% DOMAIN PARAMETERS
% domain parameters
xLim = [0 1];

% mesh parameters
base = 2;

% PDE PARAMETERS
% pde coefficients
k = 1;
r = 0;

% boundary conditions
BCtypes = 'DN';

u_N_R = -(exp(u-2) + 2);  % nonlinear Neumann BC on right edge


% BLACK BOX %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% create domain 
dom = Domain1d(xLim);

% manufacture auxiliary functions
manfun = ManufacturedFunctions1d_poisson(k,r,uTrue,BCtypes,'reaction_form','source');
auxfun = manfun.outputCoefficients();

% manufacture boundary
dom.boundary = manfun.outputBoundaryConditions(dom,BCtypes);

% overwrite appropriate boundary condition with nonlinear Neumann BC
dom.boundary.edges(2).boundaryCondition = u_N_R;
dom.boundary.edges(2).boundaryCondition_ddu = diff(u_N_R,u);


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

    % plot result
    %prob.plot;
    %pause()

end

% display results as table using matlab's table function
pList = (2:5).';
errList = cell2mat(err.');
rateList = [NaN; cell2mat(rate.')];
T = table(pList, errList, rateList);
T.Properties.VariableNames = {'p' 'L2_Error' 'Convergence_Rate'};
disp(T);
