clear all; syms x t u;
% USER INPUTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% original problem
uTrue = x + x^2;

% determine to use nonlinear BCs on left, right, or both
LEFT_nonlinearBCs = 1;      % if 1, use nonlinear BCs; if 0, use linear (manufactured) BCs
RIGHT_nonlinearBCs = 1;     % if 1, use nonlinear BCs; if 0, use linear (manufactured) BCs
constant_initial_guess = 0; % if 1, use constant initial guess u = 1; else use perturbed true solution
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
BCtypes = 'DD';

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
