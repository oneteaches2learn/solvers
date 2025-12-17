clear all; syms x t u;
% USER INPUTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% desired uTrue
%uTrue = sin(pi*x/2);

% original problem
uTrue = x + x^2;

% reversed problem
%uTrue = (1-x) + (1-x)^2;

% DOMAIN PARAMETERS
% domain parameters
xLim = [0 1];

% mesh parameters
base = 2;


% PDE PARAMETERS
% pde coefficients
k = 1;
r = u;

% boundary conditions
BCtypes = 'NN';

% original problem
U_N_L = exp(u);
dU_N_L = exp(u);
U_N_R = -(exp(u-2) + 2);
dU_N_R = -exp(u-2);

%{
% reversed problem
U_N_R = exp(u);
dU_N_R = exp(u);
U_N_L = -(exp(u-2) + 2);
dU_N_L = -exp(u-2);
%}

BCvals = {U_N_L, U_N_R};
BCvals_du = {dU_N_L, dU_N_R};

% BLACK BOX %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% create domain 
dom = Domain1d(xLim);

% manufacture auxiliary functions
manfun = ManufacturedFunctions1d_poisson(k,r,uTrue,BCtypes);
auxfun = manfun.outputCoefficients();

%{
% manually set boundary conditions (necessary for nonlinear boundary conditions)
dom.boundary = Boundary1d(xLim);
dom.boundary.bTypes = BCtypes;
dom.boundary.bConditions = BCvals;
dom.boundary.bConditions_du = BCvals_du;
%}

% manufacture boundary conditions, then replace one or both with nonlinear conditions
dom.boundary = manfun.outputBoundaryConditions(dom,BCtypes);
%dom.boundary.edges(1).boundaryCondition = U_N_L;
%dom.boundary.edges(1).boundaryCondition_ddu = dU_N_L;
dom.boundary.edges(2).boundaryCondition = U_N_R;
dom.boundary.edges(2).boundaryCondition_ddu = dU_N_R;

%{
% manufacture boundary conditions (i.e. use linear BCs)
dom.boundary = manfun.outputBoundaryConditions(dom,BCtypes);
%}

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
    uTrue_handle = matlabFunction(symfun(uTrue,x));
    prob.U = 1 * prob.U;
    %prob.U = uTrue_handle(prob.domain.mesh.nodes);
    %prob.U = prob.U + 0.1; % perturb initial guess slightly

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
