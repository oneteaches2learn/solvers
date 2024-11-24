% fem_pde_solver.m
% Discretizes the PDE:
%   u_t - d/dx (k du/dx) + r(u - u_star) = f
% on the domain Omega = (0, 1) with Dirichlet boundary conditions:
%   u(0) = u_L
%   u(1) = u_R
% This script assembles the raw stiffness matrix (K) and mass matrix (M_mat)
% using the Galerkin Finite Element Method with piecewise linear basis functions.

clear; clc;

%% Parameters
M = 20;        % Number of elements (total nodes = M + 1)
k = 1;         % Diffusion coefficient (k > 0)
r = 1;         % Reaction coefficient (r > 0)
u_star = 0;    % Given constant u^*
u_L = 0;       % Dirichlet boundary condition at x = 0
u_R = 0;       % Dirichlet boundary condition at x = 1
f_val = 0;     % Source term f (can be a function of x)

Tend = 1;      % End time for simulation

%% Spatial Discretization
x = linspace(0, 1, M+1)';      % Grid points
h = x(2) - x(1);               % Element size

% Number of interior nodes (excluding boundaries)
N = M - 1;

%% Assemble Stiffness Matrix (K) and Boundary Vector (F_boundary)
% Initialize stiffness matrix
K = sparse(N, N);       % Stiffness matrix

% Initialize boundary contribution vector
F_boundary = zeros(N, 1);

% Assemble stiffness matrix
for e = 1:M
    % Nodes associated with the current element
    n1 = e;       % Global node number
    n2 = e + 1;   % Global node number

    % Local stiffness matrix
    K_local = (k / h) * [1, -1; -1, 1];

    % Map local node indices to unknown indices
    node_indices = [n1, n2];
    unknown_indices = node_indices - 1;  % Adjust node indices to unknown indices
    local_unknowns = zeros(1, 2);

    for idx = 1:2
        if unknown_indices(idx) >= 1 && unknown_indices(idx) <= N
            local_unknowns(idx) = unknown_indices(idx);
        else
            local_unknowns(idx) = 0;  % Zero indicates a boundary node
        end
    end

    % Assemble global stiffness matrix K
    for i_local = 1:2
        i_global = local_unknowns(i_local);
        if i_global == 0
            continue;  % Skip boundary nodes
        end
        for j_local = 1:2
            j_global = local_unknowns(j_local);
            if j_global == 0
                continue;  % Skip boundary nodes
            end
            K(i_global, j_global) = K(i_global, j_global) + K_local(i_local, j_local);
        end
    end
end


%% Assemble Mass Matrix (M_mat)
% Initialize mass matrix
M_mat = sparse(N, N);   % Mass matrix

% Assemble mass matrix
for e = 1:M
    % Nodes associated with the current element
    n1 = e;
    n2 = e + 1;

    % Local mass matrix
    M_local = (h / 6) * [2, 1; 1, 2];

    % Map local node indices to unknown indices
    node_indices = [n1, n2];
    unknown_indices = node_indices - 1;  % Adjust node indices to unknown indices
    local_unknowns = zeros(1, 2);

    for idx = 1:2
        if unknown_indices(idx) >= 1 && unknown_indices(idx) <= N
            local_unknowns(idx) = unknown_indices(idx);
        else
            local_unknowns(idx) = 0;  % Zero indicates a boundary node
        end
    end

    % Assemble global mass matrix M_mat
    for i_local = 1:2
        i_global = local_unknowns(i_local);
        if i_global == 0
            continue;  % Skip boundary nodes
        end
        for j_local = 1:2
            j_global = local_unknowns(j_local);
            if j_global == 0
                continue;  % Skip boundary nodes
            end
            M_mat(i_global, j_global) = M_mat(i_global, j_global) + M_local(i_local, j_local);
        end
    end
end

%% Source Term Vector (F)
% Assuming f is constant; modify if f is a function of x
F = f_val * ones(N, 1);

%% Initial Condition
% Initial condition at interior nodes (excluding boundaries)
U0 = zeros(N, 1);
% If you have an initial condition u0(x), evaluate it at interior nodes:
% U0 = u0(x(2:end-1));

%% ODE System Setup
% Mass matrix for the reaction term
M_react = r * M_mat;

% Total mass matrix (diffusion and reaction combined)
M_total = M_mat;

% Ensure mass matrix is full
M_total = full(M_total);

% Combine stiffness and reaction matrices
A = K + M_react;

% Convert matrices to full to avoid dimension issues
A = full(A);
F_total = F_boundary + F + M_react * u_star;

% Right-hand side function
RHS = @(t, U) F_total - A * U;
A

%% Solve ODE System using ode15s
% Define options with mass matrix
%opts = odeset('Mass', M_total, 'MStateDependence', 'none', 'MassSingular', 'no');

% Time span
tspan = [0, Tend];

% Solve the system
%[t_sol, U_sol] = ode45(@(t, U) RHS(t, U), tspan, U0, opts);
[t_sol, U_sol] = ode45(@(t, U) RHS(t, U), tspan, U0);

return
%% Post-processing and Visualization
% Include boundary values
U_full = [u_L * ones(length(t_sol), 1), U_sol, u_R * ones(length(t_sol), 1)];

% Plot solution at final time
figure;
plot(x, U_full(end, :), 'b-o', 'LineWidth', 1.5);
xlabel('Spatial Variable x');
ylabel('Solution u(x, t_{end})');
title('Solution at Final Time');
grid on;

% Surface plot over time
figure;
surf(x, t_sol, U_full, 'EdgeColor', 'none');
xlabel('Spatial Variable x');
ylabel('Time t');
zlabel('Solution u(x, t)');
title('Solution Surface over Space and Time');
colorbar;
view(130, 30);