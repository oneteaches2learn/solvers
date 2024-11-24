% galerkin_1d.m
% Galerkin Finite Element Solver for the Poisson Equation:
%   -k * d^2u/dx^2 = f(x) on Omega = (0,1)
% with homogeneous Dirichlet boundary conditions:
%   u(0) = 0 and u(1) = 0

clear; clc; close all;

%% 1. Parameters and Domain Discretization

% Physical Parameters
k = 1;                 % Diffusion coefficient (k > 0)
f = @(x) sin(pi*x);    % Source term f(x)

% Numerical Parameters
M = 20;                % Number of finite elements (total nodes = M + 1)
L = 1;                 % Length of the domain
x = linspace(0, L, M+1)';  % Grid points (column vector)
h = x(2) - x(1);           % Element size

% Number of Interior Nodes (excluding boundaries)
N = M - 1;

%% 2. Assembly of Stiffness Matrix (K) and Load Vector (F)

% Initialize global stiffness matrix and load vector
K = sparse(N, N);   % Stiffness matrix
F = zeros(N, 1);    % Load vector

% Loop over each finite element to assemble K and F
for e = 1:M
    % Nodes associated with the current element
    n1 = e;     % Global node number (1 to M)
    n2 = e + 1; % Global node number (2 to M+1)
    
    % Local stiffness matrix for the element
    K_local = (k / h) * [1, -1; -1, 1];
    
    % Local load vector for the element (using Trapezoidal Rule)
    F_local = (h / 6) * [2*f(x(n1)) + f(x(n2)); f(x(n1)) + 2*f(x(n2))];
    % Alternatively, for higher accuracy, use Simpson's Rule or Gaussian Quadrature
    
    % Map local nodes to global unknowns
    % Interior nodes are 2 to M, corresponding to unknowns 1 to N
    % Nodes 1 and M+1 are boundary nodes with known values (homogeneous Dirichlet)
    global_nodes = [n1, n2];
    global_unknowns = global_nodes - 1;
    
    % Assemble the global stiffness matrix
    for i = 1:2
        for j = 1:2
            % Only assemble if both nodes are interior
            if (global_unknowns(i) >= 1) && (global_unknowns(i) <= N) && ...
               (global_unknowns(j) >= 1) && (global_unknowns(j) <= N)
                K(global_unknowns(i), global_unknowns(j)) = ...
                    K(global_unknowns(i), global_unknowns(j)) + K_local(i, j);
            end
        end
    end
    
    % Assemble the global load vector
    for i = 1:2
        if (global_unknowns(i) >= 1) && (global_unknowns(i) <= N)
            F(global_unknowns(i)) = F(global_unknowns(i)) + F_local(i);
        end
    end
end

%% 3. Application of Homogeneous Dirichlet Boundary Conditions

% Since boundary conditions are homogeneous (u(0)=0 and u(1)=0),
% and we have excluded boundary nodes from the unknowns, no modification
% to K and F is necessary. If boundary conditions were non-homogeneous,
% adjustments to F would be required here.

%% 4. Solution of the Linear System

% Solve for interior node values
u_interior = K \ F;

% Construct the full solution including boundary nodes
u = [0; u_interior; 0];

%% 5. Post-processing and Visualization

% Analytical solution (optional, for verification)
% For f(x) = sin(pi x), the analytical solution is:
% u(x) = (1 / pi^2) * sin(pi x)
u_exact = (1 / pi^2) * sin(pi*x);

% Plot Numerical and Analytical Solutions
figure;
plot(x, u, 'bo-', 'LineWidth', 1.5, 'MarkerFaceColor', 'b');
hold on;
plot(x, u_exact, 'r--', 'LineWidth', 2);
xlabel('Spatial Coordinate x');
ylabel('Solution u(x)');
title('Finite Element Solution vs. Analytical Solution');
legend('Numerical FEM Solution', 'Analytical Solution');
grid on;
hold off;

% Plot Absolute Error
error = abs(u - u_exact);
figure;
plot(x, error, 'k-o', 'LineWidth', 1.5, 'MarkerFaceColor', 'k');
xlabel('Spatial Coordinate x');
ylabel('|u_{numerical} - u_{exact}|');
title('Absolute Error Between Numerical and Analytical Solutions');
grid on;

%% 6. (Optional) Display Assembled Matrices

% Uncomment below lines to display the stiffness matrix and load vector
% disp('Stiffness Matrix K:');
% full(K)

% disp('Load Vector F:');
% F
