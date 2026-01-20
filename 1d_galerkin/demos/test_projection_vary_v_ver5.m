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
p = 2;


% PDE PARAMETERS
% pde coefficients
lambda = 10;
k = 1e-4;
r = lambda * u;
dr_du = lambda;

% boundary conditions
BCtypes = 'NN';

% original problem
%U_N = - v^3 - u^3;
%dU_N = diff(U_N,u);

beta = 5;
v0   = 1;
U_N  = -(tanh(v/v0)) * tanh(beta*u);   % bounded in both u and v
dU_N = diff(U_N,u);

BCvals = {U_N, U_N};
BCvals_du = {dU_N, dU_N};

% BLACK BOX %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% create domain 
dom = Domain1d(xLim);

% make waitbar
%wait = waitbar(0,'Solving projection problems...');

% initialize solution storage
ind = 1;
exponents_pos = 10.^(-5:1:2);
exponents_neg = -fliplr(exponents_pos);
exponents = [exponents_neg 0 exponents_pos];

N = length(exponents);

hold on
for i = exponents
%for i = 10
%for i = -2:.1:0

    % set u and v
    %u_fixed = x + x^2 + 1;
    %u_fixed = 1 + x + x^2;
    %u_fixed = sin(2 * pi * x);
    %u_fixed = x^2;
    u_fixed = exp(5*x);
    v_fixed = i;

    % manufacture auxiliary functions
    auxfun = AuxFun1d_projection(k,0);
    auxfun.r = matlabFunction(symfun(r,x));
    auxfun.dr_du = matlabFunction(symfun(dr_du,x));
    auxfun.u_fixed = matlabFunction(symfun(u_fixed,x));
    auxfun.v_fixed = matlabFunction(symfun(v_fixed,x));
    auxfun.lambda = lambda;

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

    % store solution
    sol(ind,:) = prob.solution';
    ind = ind + 1;

    % --- Quick L2 error (nodal + trapz) ---
    xnodes = prob.domain.mesh.nodes(:);                 % column
    u_vals = auxfun.u_fixed(xnodes);                    % u(x) at nodes (column)
    w_vals = prob.solution(:);                          % W_h at nodes (column)

    L2_err(ind) = sqrt(trapz(xnodes, (w_vals - u_vals).^2));

    fprintf(' ind: %d, v = %.5e, W_max = %.3f, ||W-u||_L2 ~= %.3e\n', ...
        ind, i, max(abs(w_vals)), L2_err(ind));

    % update waitbar in a simple way
    %waitbar(ind/N,wait);

    %fprintf(' ind: %d, v = %.5f, W_max = %.3f \n',ind,i,max(abs(prob.solution)));
    plot(prob.domain.mesh.nodes,prob.solution,'LineWidth',1);
    pause()

    %{
    % compute exact solution
    u_exact_denom = 24 * (-lambda^2 + 8*lambda*v_fixed - 48*lambda + 96*v_fixed);
    u_exact(1,1) = lambda * (lambda + 48);
    u_exact(2,1) = -5 * lambda^2 + 36 * lambda * v_fixed - 240 * lambda + 576 * v_fixed;
    u_exact(3,1) = -23 * lambda^2 + 192 * lambda * v_fixed - 1104 * lambda + 2304 * v_fixed;
    u_exact = u_exact / u_exact_denom;
    %}

    %{
    % compute numerical solution
    K = [2 -2 0; -2 4 -2; 0 -2 2];
    M = [1/6 1/12 0; 1/12 1/3 1/12; 0 1/12 1/6];
    B = [1 0 0; 0 0 0; 0 0 1];
    b = [-1/2; -1; 3/2] + lambda * [1/96; 7/48; 17/96] - v_fixed * [0; 0 ;1];

    % assemble
    A = K + lambda * M - v_fixed * B;

    % solve
    w = A \ b;
    %}

    %{
    % compute error
    max(prob.solution - u_exact)
    %u_exact - w
    pause()
    %}

end
hold off

% close waitbar
%close(wait);