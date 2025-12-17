clear all; syms x t;
% USER INPUTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% desired uTrue
uTrue = sin(pi*x / 4);


% DOMAIN PARAMETERS
% domain parameters
xLim = [0 1];

% mesh parameters
base = 2;


% PDE PARAMETERS
% pde coefficients
k = 2 + sin(pi*x);
r = 1 + x^2;

% boundary conditions
BCtypes = 'NN';


% BLACK BOX %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% create domain 
dom = Domain1d(xLim);

% manufacture auxiliary functions
manfun = ManufacturedFunctions1d_poisson(k,r,uTrue,BCtypes);
auxfun = manfun.outputCoefficients();

% manufacture boundary conditions
BCvals = manfun.outputBoundaryConditions();
dom.boundary = Boundary1d(xLim,BCtypes,BCvals);

% loop over mesh refinements
err = {};
i = 0;
for p = 2:5

    % update counter
    i = i + 1;

    % update mesh
    dom.mesh = Mesh1d(xLim,base,p);

    % create galerkin solver
    prob = GalerkinSolver1d_poisson(dom,auxfun);
    prob = prob.solve;

    % compute and display errors
    err{i} = dom.L2err_twoPointQuadrature(prob.solution,uTrue);

    % compute difference in errors
    if i > 1
        rate{i-1} = log(err{i-1}/err{i})/log(base);
    end

end

% display results as table using matlab's table function
pList = (2:5).';
errList = cell2mat(err.');
rateList = [NaN; cell2mat(rate.')];
T = table(pList, errList, rateList);
T.Properties.VariableNames = {'p' 'L2_Error' 'Convergence_Rate'};
disp(T);
