% galerkin1d_semidiscrete_coupled_dirichlet.m
% Solves the time-dependent 1D FEM problem using GalerkinSolver1d and ode45
% Includes an animation of the solution over time

clear; clc; close all;

%% USER INPUTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PDE parameters
k  = 1;  % Diffusion coefficient (k > 0)
r  = 1;  % Reaction coefficient (r > 0)
f  = 1;
u0 = 0;

% ODE parameters
s  = 1;
g  = 0;
v0 = 0;

% Domain Parameters
M = 100;       % Number of finite elements
L = 1;       % Length of the domain
T_end = 4;   % End time for ODE solver
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
U0 = u0_full(2:end-1);  % Interior nodes

% Compute total initial condition
W0 = [U0; v0];

% Define the ODE function
% PDE contribution to system
K = elliptic.tensors.K;
R = elliptic.tensors.R;
F = cofs.f(domain.x(2:end-1));
Mat = K + R;

% ODE contribution to system
rV = r * ones(size(K,1),1);
SU = s * domain.h * ones(1,size(K,2));
V = s;
G = g;

% Monolithic system
Fbar = [F; G];
Mbar = [Mat, -rV; -SU, V];

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
U_full = zeros(size(U,1), size(U,2)+2);
U_full(:,2:end-1) = U;

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