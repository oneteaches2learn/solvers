% galerkin1d_semidiscrete_coupled.m
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
g  = 1;
v0 = 0;

% Domain Parameters
M = 10;      % Number of finite elements
L = 1;       % Length of the domain
T_end = 1;   % End time for ODE solver

%% BLACK BOX %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Store coefficients in structure
cofs.k = k;                 % Diffusion coefficient
cofs.r = r;                 % Reaction coefficient
cofs.f = @(x) 1;    % Source term f(x)

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
F = elliptic.vectors.F;
Mat = K + R;

% ODE contribution to system
rV = r * ones(size(K,1),1);
SU = s * ones(1,size(K,2));
V = s;
G = 0;

% Monolithic system
Fbar = [F; G];
Mbar = [Mat, rV; SU, V];

% final ode system
odefun = @(t, W) Fbar - Mbar * W;

%% SOLVE ODE USING ode45 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Time span
tspan = [0, T_end];

% Solve the ODE
[t, W] = ode45(odefun, tspan, W0);


%% POST-PROCESSING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Separate U and V
U = W(:,1:end-1);
V = W(:,end);

% Add back homogeneous Dirichlet BCs
U_full = zeros(size(U,1), size(U,2)+2);
U_full(:,2:end-1) = U;
size(U_full)
size(domain.x)

% Compute y-limits for Subplot 1 (All columns except last)
min_u = min(U_full(:,:), [], 'all');
max_u = max(U_full(:,:), [], 'all');
range_u = max_u - min_u;
ymin_u = min_u - 0.1 * range_u;
ymax_u = max_u + 0.1 * range_u;

% Compute y-limits for Subplot 2 (Last column over time)
min_v = min(V);
max_v = max(V);
range_v = max_v - min_v;
ymin_v = min_v - 0.1 * range_v;
ymax_v = max_v + 0.1 * range_v;

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
    
    pause(0.01);  % Adjust pause duration as needed
end
hold off;

