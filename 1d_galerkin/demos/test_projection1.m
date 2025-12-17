% test_projection1.m
% Loops over different values of v_fixed and approximates the solution to the
% projection problem, then plots these solutions on the same graph. 

clear all; syms x t u v;
% USER INPUTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% original problem


% DOMAIN PARAMETERS
% domain parameters
xLim = [0 1];

% mesh parameters
base = 2;
p = 5;


% PDE PARAMETERS
% pde coefficients
lambda = 100;
k = 1;
r = lambda * u;
dr_du = lambda;

% boundary conditions
BCtypes = 'NN';

% original problem
U_N_L = u^3 * v + v^3;
U_N_R = u^3 * v + v^3;
dU_N_L = diff(U_N_L,u);
dU_N_R = diff(U_N_R,u);

BCvals = {U_N_L, U_N_R};
BCvals_du = {dU_N_L, dU_N_R};

% BLACK BOX %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% create domain 
dom = Domain1d(xLim);

% make waitbar
wait = waitbar(0,'Solving projection problems...');

hold on
for i = -500000:100000:500000

    % set u and v
    u_fixed = x + x^2 + 1;
    v_fixed = i;

    % update waitbar in a simple way
    %waitbar((i+50000)/100000,wait);

    % manufacture auxiliary functions
    auxfun = AuxFun1d_projection(k,0);
    auxfun.r = matlabFunction(symfun(r,x));
    auxfun.dr_du = matlabFunction(symfun(dr_du,x));
    auxfun.u_fixed = matlabFunction(symfun(u_fixed,x));
    auxfun.v_fixed = matlabFunction(symfun(v_fixed,x));

    % add boundary conditions to domain
    dom.boundary = Boundary1d(xLim);
    dom.boundary.bTypes = BCtypes;
    dom.boundary.bConditions = BCvals;
    dom.boundary.bConditions_du = BCvals_du;


    % mesh the domain
    dom.mesh = Mesh1d(xLim,base,p);

    % create the projection problem and solve
    prob = NewtonGalerkinSolver1d_projection(dom,auxfun);


    % solve the problem with the updated v_fixed
    prob = prob.solve;

    plot(prob.domain.mesh.nodes,prob.solution,'LineWidth',1);

end
hold off

% close waitbar
close(wait);