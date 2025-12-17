% USER INPUTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% original problem
clear all;
hold on
for i = 1:5
syms x t u v;
u_fixed = 1 + x^2;
v_fixed = -1 - .1*i;


% DOMAIN PARAMETERS
% domain parameters
xLim = [0 1];

% mesh parameters
base = 2;
p = 5;


% PDE PARAMETERS
% pde coefficients
lambda = 1;
k = 1;
r = lambda * u;
dr_du = lambda;

% boundary conditions
BCtypes = 'NN';

% original problem
U_N_L = u^3 + u*v;
U_N_R = u^3 + u*v;
dU_N_L = diff(U_N_L,u);
dU_N_R = diff(U_N_R,u);

BCvals = {U_N_L, U_N_R};
BCvals_du = {dU_N_L, dU_N_R};

% BLACK BOX %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% create domain 
dom = Domain1d(xLim);

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

% overwrite b_vol with exact value computed by hand
%b_vol_exact = [-59/48 - v_fixed; 19/12; -59/48 - v_fixed];
b_vol_exact = [-119/96 - v_fixed; -17/48; -583/96 - 2*v_fixed];
%prob.vectors.b_vol = b_vol_exact;

% solve
prob = prob.solve;

% plot the solution
plot(prob.domain.mesh.nodes,prob.solution,'LineWidth',2);
pause()

clear all
end


%{
tensors = prob.tensors;
b_stiff = [-1/2; -1; 3/2];
b_mass  = [1/96; 7/48; 17/96];
b_neu   = [0; 0; 1];
b_vol   = b_stiff + lambda * b_mass - b_neu;
sol_exact = (tensors.A + tensors.M_r - tensors.M_dneu) \ b_vol;
%}







%{
% print results
b_vol_exact = b_vol
b_vol_computed = prob.vectors.b_vol
sol_exact
prob.solution
%}


%{
hold on
for i = -1000:100:1000

    prob.coefficients.v_fixed = matlabFunction(symfun(i-5,x));
    %prob.coefficients.lambda = i;
    prob = prob.solve;

    % plot the solution
    plot(prob.domain.mesh.nodes,prob.solution)

end
hold off
%}
