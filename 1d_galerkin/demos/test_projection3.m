% test_projection1.m
% Loops over different values of v_fixed and approximates the solution to the
% projection problem, then plots these solutions on the same graph. 

clear all; syms x t u v;
% USER INPUTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% original problem
clf;
hold on
for i = [-10,0.1,10]
u_fixed = 1 + x - x^2;
v_fixed = i;


% DOMAIN PARAMETERS
% domain parameters
xLim = [0 1];

% mesh parameters
base = 2;
p = 1;


% PDE PARAMETERS
% pde coefficients
lambda = 1;
k = 1;
r = lambda * u;
dr_du = lambda;

% boundary conditions
BCtypes = 'NN';

% original problem
U_N_L = u^3 * v;
U_N_R = u^3 * v;
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

% overwrite b_vol
b_vol_exact = [-1/2 + 9/32 - v_fixed; 1 + 29/48; -1/2 + 9/32 - v_fixed];
prob.vectors.b_vol = b_vol_exact;

% solve and plot
prob = prob.solve;
plot(prob.domain.mesh.nodes,prob.solution,'LineWidth',2);


end
hold off

f = gcf
g = f.Children;
g.Children(1).LineWidth = 8;
g.Children(1).LineStyle = '--';
g.Children(1).Color = [1 0 0];

g.Children(2).LineWidth = 6;
g.Children(2).Color = [0 0 1];

g.Children(3).LineWidth = 6;
g.Children(3).Color = [0 0 0];

g.FontSize = 24;

legend('v_{fixed} = -10','v_{fixed} = 0.1','v_{fixed} = 10');
