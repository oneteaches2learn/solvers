% galerkin1d_semidiscrete_coupled_dirichlet.m
% Solves the time-dependent 1D FEM problem using GalerkinSolver1d and ode45
% Includes an animation of the solution over time

clear; clc; close all;
%% USER INPUTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PDE parameters
k  = 1;  % Diffusion coefficient (k > 0)
r  = 1;  % Reaction coefficient (r > 0)
f  = 0;
u0 = 2;

% ODE parameters
s  = 1;
g  = 0;
v0 = 2;

% Domain Parameters
M = 10;       % Number of finite elements
L = 1;       % Length of the domain
T_end = 10;   % End time for ODE solver
Nt = 100;    % Number of time steps

%% BLACK BOX %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Store coefficients in structure
cofs.k = k;                 % Diffusion coefficient
cofs.r = r;                 % Reaction coefficient
cofs.f = @(x) 0*x + f;    % Source term f(x)

% Store domain information in structure
domain.M = M;                                    % Number of finite elements
domain.L = L;                                    % Length of the domain
domain.x = linspace(0, domain.L, domain.M + 1)'; % Grid points (column vector)
domain.h = domain.x(2) - domain.x(1);            % Element size
domain.N = domain.M - 1;                         % Number of Interior Nodes (excluding boundaries)

% Create elliptic problem to generate tensors
elliptic = GalerkinSolver1d(domain, cofs);

%% DEFINE ODE SYSTEM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute PDE initial condition
u0 = @(x) 0*x + u0;
u0_full = elliptic.evaluateFunction(u0);
U0 = u0_full;

% Compute total initial condition
W0 = [U0; v0];

% Define the ODE function
% PDE contribution to system
K = assembleStiffnessMatrix(domain.M, domain.N, domain.h, k);
R = assembleMassMatrix(domain.M, domain.N, domain.h, r);
F = cofs.f(domain.x);
Mat = K + R;

% ODE contribution to system
rV = r * domain.h * ones(size(K,1),1);
SU = s * domain.h * ones(1,size(K,2));
SU(1) = SU(1)/2;
SU(end) = SU(end)/2;
V = s;
G = g;

% neumann boundary conditions
F_neu = generateNeumannBoundaryConditions(F, 10, 10, domain.h);
%F = F + F_neu;
rV = rV - F_neu;


% Monolithic system
Fbar = [F; G];
Mbar = [Mat, -rV; -SU, V];

fprintf('Max and min eigenvalues of system:\n\n')
Eig = eig(full(-Mbar));
max_eig = max(Eig)
min_eig = min(Eig)

% final ode system
odefun = @(t, W) Fbar - Mbar * W;

%% SOLVE ODE USING ode45 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Time span
tspan = linspace(0,T_end,Nt);

% Solve the ODE
[t, W] = ode45(odefun, tspan, W0);

%% POST-PROCESSING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Separate U and V
U = W(:,1:end-1);
V = W(:,end);

% Add back homogeneous Dirichlet BCs
U_full = U;

% compute U_avg
U_avg = domain.h * (U_full(:,1)/2 + U_full(:,end)/2 + sum(U_full(:,2:end-1),2));

% Compute y-limits for Subplot 1 (All columns except last)
min_u = min(U_full(:,:), [], 'all');
max_u = max(U_full(:,:), [], 'all');
range_u = max_u - min_u;
ymin_u = min_u - 0.1 * range_u;
ymax_u = max_u + 0.1 * range_u;
if ymin_u == ymax_u
    ymin_u = ymin_u - 0.1;
    ymax_u = ymax_u + 0.1;
end

% Compute y-limits for Subplot 2 (Last column over time)
min_v = min(V);
max_v = max(V);
range_v = max_v - min_v;
ymin_v = min_v - 0.1 * range_v;
ymax_v = max_v + 0.1 * range_v;
if ymin_v == ymax_v
    ymin_v = ymin_v - 0.1;
    ymax_v = ymax_v + 0.1;
end

% Compute limits for U_avg
min_uAvg = min(U_avg);
max_uAvg = max(U_avg);
range_uAvg = max_uAvg - min_uAvg;
ymin_uAvg = min_uAvg - 0.1 * range_uAvg;
ymax_uAvg = max_uAvg + 0.1 * range_uAvg;
if ymin_uAvg == ymax_uAvg
    ymin_uAvg = ymin_uAvg - 0.1;
    ymax_uAvg = ymax_uAvg + 0.1;
end

%{
%% CREATE ANIMATION WITH SUBPLOTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
figure('Position', [100, 100, 800, 600]);

% Subplot 1: FEM Solution Evolution (All Columns Except Last)
subplot(1,2,1);
h_anim = plot(domain.x(1:end), U_full(1,:), 'b-', 'LineWidth', 1.5);
hold on;
xlabel('Spatial Coordinate x');
ylabel('Solution u(x, t)');
title('FEM Solution Evolution (Excl. Last Column)');
legend('FEM Solution');
grid on;
ylim([ymin_u, ymax_u]); % Set y-limits
hold off;

% Subplot 2: Last Column Over Time
subplot(1,2,2);
h_v = plot(t(1), V(1), 'k-', 'LineWidth', 1.5);
xlabel('Time t');
ylabel('v(t)');
title('Coupled Variable v(t) Over Time');
legend('v(t)');
grid on;
hold on;
ylim([ymin_v, ymax_v]); % Set y-limits
xlim([0, T_end]);       % Set x-limits

% Animate the plots
pause();
for i = 1:length(t)
    % Update FEM Solution in Subplot 1
    subplot(1,2,1);
    set(h_anim, 'YData', U_full(i,:));
    drawnow;
    
    % Update v(t) in Subplot 2
    subplot(1,2,2);
    set(h_v, 'XData', t(1:i));
    set(h_v, 'YData', V(1:i));
    drawnow;
    
    pause(1/length(t));  % Adjust pause duration as needed
end
hold off;
%}


%% CREATE PHASE PLOT AND SUPERIMPOSED PLOTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure('Position', [150, 150, 1200, 600]);  % [left, bottom, width, height]

% --- Subplot 1: Phase Plot (U_avg vs V) ---
subplot(1,2,1);
hold on
h_phase = plot(U_avg, V, '-', 'Color','b','LineWidth',1.5);
plot(U_avg(1), V(1), 'ko', 'MarkerFaceColor','g');
xlabel('U_{avg}(t)');
ylabel('V(t)');
title('Phase Plot: U_{avg} vs V');
legend('Trajectory','Initial Condition');
grid on;
ylim([ymin_v, ymax_v]);
xlim([ymin_uAvg, ymax_uAvg]);

% --- Subplot 2: U_avg and V Superimposed Over Time ---
subplot(1,2,2);
h_U_avg = plot(t, U_avg, 'b-', 'LineWidth', 1.5);
hold on;
h_V = plot(t, V, 'r-', 'LineWidth', 1.5);
xlabel('Time t');
ylabel('Values');
title('U_{avg} and V over Time');
legend('U_{avg}(t)', 'V(t)');
grid on;
xlim([0, T_end]);
temp_min = min([ymin_uAvg, ymin_v]);
temp_max = max([ymax_uAvg, ymax_v]); 
temp_range = temp_max - temp_min;
ylim([temp_min - 0.1*temp_range, temp_max + 0.1*temp_range]);

%{
% --- Animate the Phase Plot and Superimposed Plots ---
pause()
for i = 1:length(t)
    % Update Phase Plot in Subplot 1
    subplot(1,2,1);
    set(h_phase, 'XData', U_avg(1:i));
    set(h_phase, 'YData', V(1:i));
    drawnow;
    
    % Update Superimposed Plots in Subplot 2
    subplot(1,2,2);
    set(h_U_avg, 'XData', t(1:i));
    set(h_U_avg, 'YData', U_avg(1:i));
    set(h_V, 'XData', t(1:i));
    set(h_V, 'YData', V(1:i));
    drawnow;
    
    %pause(1/length(t));  % Adjust pause duration as needed
end
hold off;
%}


%{
%% CREATE PHASE ANIMATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure('Position', [150, 150, 1200, 600]);  % [left, bottom, width, height]

% --- Subplot 1: Phase Plot (U_avg vs V) ---
subplot(1,2,1);
hold on
h_phase = plot(U_avg(1), V(1), '-', 'Color','b','LineWidth',1.5);
plot(U_avg(1), V(1), 'ko', 'MarkerFaceColor','g');
xlabel('U_{avg}(t)');
ylabel('V(t)');
title('Phase Plot: U_{avg} vs V');
legend('Trajectory','Initial Condition');
grid on;
ylim([ymin_v, ymax_v]);
xlim([ymin_uAvg, ymax_uAvg]);

% --- Subplot 2: U_avg and V Superimposed Over Time ---
subplot(1,2,2);
h_U_avg = plot(t(1), U_avg(1), 'b-', 'LineWidth', 1.5);
hold on;
h_V = plot(t(1), V(1), 'r-', 'LineWidth', 1.5);
xlabel('Time t');
ylabel('Values');
title('U_{avg} and V over Time');
legend('U_{avg}(t)', 'V(t)');
grid on;
xlim([0, T_end]);
temp_min = min([ymin_uAvg, ymin_v]);
temp_max = max([ymax_uAvg, ymax_v]); 
temp_range = temp_max - temp_min;
ylim([temp_min - 0.1*temp_range, temp_max + 0.1*temp_range]);

% --- Animate the Phase Plot and Superimposed Plots ---
pause()
for i = 1:length(t)
    % Update Phase Plot in Subplot 1
    subplot(1,2,1);
    set(h_phase, 'XData', U_avg(1:i));
    set(h_phase, 'YData', V(1:i));
    drawnow;
    
    % Update Superimposed Plots in Subplot 2
    subplot(1,2,2);
    set(h_U_avg, 'XData', t(1:i));
    set(h_U_avg, 'YData', U_avg(1:i));
    set(h_V, 'XData', t(1:i));
    set(h_V, 'YData', V(1:i));
    drawnow;
    
    %pause(1/length(t));  % Adjust pause duration as needed
end
hold off;
%}


%% ASSEMBLE STIFFNESS MATRIX
function K = assembleStiffnessMatrix(M, N, h, k)
    % Assembles the global stiffness matrix K for homogeneous Neumann boundary conditions
    
    % Total number of nodes (including boundary nodes)
    totalNodes = N + 2;  % N interior nodes + 2 boundary nodes
    
    % Initialize global stiffness matrix
    K = sparse(totalNodes, totalNodes);
    
    % Local stiffness matrix for each element
    K_local = (k / h) * [1, -1; -1, 1];
    
    % Assemble global stiffness matrix
    for e = 1:M
        % Nodes associated with the current element
        n1 = e;       % Global node number (1 to M)
        n2 = e + 1;   % Global node number (2 to M+1)
        
        % Map local nodes to global nodes (no adjustment needed)
        global_nodes = [n1, n2];
        
        % Assemble the stiffness matrix
        for i = 1:2
            for j = 1:2
                K(global_nodes(i), global_nodes(j)) = ...
                    K(global_nodes(i), global_nodes(j)) + K_local(i, j);
            end
        end
    end
end

%% ASSEMBLE MASS MATRIX
function R = assembleMassMatrix(M,N,h,c)
    % Assembles the global mass matrix R scaled by coefficient c
    
    % Total number of nodes (including boundary nodes)
    totalNodes = N + 2;  % N interior nodes + 2 boundary nodes

    % Initialize global mass matrix
    R = sparse(totalNodes, totalNodes);
    
    % Local mass matrix for each element
    R_local = c * (h / 6) * [2, 1; 1, 2];
    
    % Assemble global mass matrix
    for e = 1:M
        % Nodes associated with the current element
        n1 = e;         % Global node number (1 to M)
        n2 = e + 1;     % Global node number (2 to M+1)
        
        % Map local nodes to global nodes (no adjustment needed)
        global_nodes = [n1, n2];
        
        % Assemble the mass matrix
        for i = 1:2
            for j = 1:2
                % Add contributions to the global mass matrix
                R(global_nodes(i), global_nodes(j)) = ...
                    R(global_nodes(i), global_nodes(j)) + R_local(i, j);
            end
        end
    end
end


function F_neu = generateNeumannBoundaryConditions(F, q0, qL, h)
    % GENERATENEUMANNBOUNDARYCONDITIONS Generates a F_neu to account for
    % inhomogeneous Neumann boundary conditions at the boundaries x = 0 and x = L.
    %
    % Syntax:
    %   F_modified = applyNeumannBoundaryConditions(F, q0, qL)
    %
    % Inputs:
    %   F  - Existing load vector (column vector of size [totalNodes x 1])
    %   q0 - Neumann boundary condition at x = 0 (flux or derivative value)
    %   qL - Neumann boundary condition at x = L (flux or derivative value)
    %
    % Output:
    %   F_modified - Modified load vector including Neumann boundary contributions
    %
    % Notes:
    %   - The function assumes that positive flux is defined outward from the domain.
    %   - The Neumann conditions q0 and qL should be provided based on the problem's
    %     physical context and sign conventions.

    F_neu = sparse(size(F,1), 1);  % Initialize modified load vector

    % Add Neumann boundary contributions to the load vector
    %F_neu(1)   = q0 * h / 2;  % Contribution at x = 0
    %F_neu(end) = qL * h / 2;  % Contribution at x = L
    F_neu(1)   = q0 / 2;  % Contribution at x = 0
    F_neu(end) = qL / 2;  % Contribution at x = L
end