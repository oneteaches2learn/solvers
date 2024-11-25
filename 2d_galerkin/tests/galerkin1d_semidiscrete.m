% galerkin1d_semidiscrete.m
% Solves the time-dependent 1D FEM problem using GalerkinSolver1d and ode45
% Includes an animation of the solution over time

clear; clc; close all;

%% USER INPUTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Desired true solution
syms x
u_true_sym(x) = sin(pi * x);  % User can modify this as desired

% Physical Parameters
k = 1;  % Diffusion coefficient (k > 0)
r = 1;  % Reaction coefficient (r > 0)

% Domain Parameters
M = 10;      % Number of finite elements
L = 1;       % Length of the domain
T_end = 1;   % End time for ODE solver

%% BLACK BOX %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute source term
f_sym = -k * diff(u_true_sym, x, 2) + r * u_true_sym;

% Convert to function handles
uTrue_func = matlabFunction(u_true_sym, 'Vars', x);
f_func = matlabFunction(f_sym, 'Vars', x);

% Store coefficients in structure
cofs.k = k;                 % Diffusion coefficient
cofs.r = r;                 % Reaction coefficient
cofs.f = @(x) f_func(x);    % Source term f(x)

fprintf('Source Function f(x) = %s\n', char(f_sym));

% Store domain information in structure
domain.M = M;                                    % Number of finite elements
domain.L = L;                                    % Length of the domain
domain.x = linspace(0, domain.L, domain.M + 1)'; % Grid points (column vector)
domain.h = domain.x(2) - domain.x(1);            % Element size
domain.N = domain.M - 1;                         % Number of Interior Nodes (excluding boundaries)

% Create elliptic problem to generate tensors
elliptic = GalerkinSolver1d(domain, cofs);

%% DEFINE ODE SYSTEM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute initial condition using u_true_sym
u0_full = elliptic.evaluateFunction(u_true_sym(x));
U0 = u0_full(2:end-1);  % Interior nodes

% Define the ODE function
K = elliptic.tensors.K;
R = elliptic.tensors.R;
F = elliptic.vectors.F;
odefun = @(t, U) F - (K + R) * U;

size(K)
%% SOLVE ODE USING ode45 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Time span
tspan = [0, T_end];

% Solve the ODE
[t, U] = ode45(odefun, tspan, U0);

%% POST-PROCESSING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Assemble full solution including boundary nodes for each time step
num_time_steps = length(t);
U_full = zeros(domain.M + 1, num_time_steps);
for i = 1:num_time_steps
    U_full(:, i) = [0; U(i, :)'; 0];
end

%% CREATE ANIMATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
h_anim = plot(domain.x, U_full(:,1), 'b-', 'LineWidth', 1.5);
hold on;
h_true = plot(domain.x, uTrue_func(domain.x), 'r--', 'LineWidth', 1.5);
xlabel('Spatial Coordinate x');
ylabel('Solution u(x, t)');
title('FEM Solution Animation');
legend('FEM Solution', 'True Solution');
grid on;
hold on;

pause()
for i = 1:num_time_steps
    set(h_anim, 'YData', U_full(:,i));
    drawnow;
    pause(0.01);  % Adjust pause duration as needed
end
hold off;

%{
%% PLOT SOLUTION AT FINAL TIME %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
plot(domain.x, U_full(:, end), 'b-', 'LineWidth', 1.5);
hold on;
plot(domain.x, uTrue_func(domain.x), 'r--', 'LineWidth', 1.5);
xlabel('Spatial Coordinate x');
ylabel('Solution u(x, T_{end})');
title('FEM Solution at Final Time vs. True Solution');
legend('FEM Solution', 'True Solution');
grid on;
hold off;
%}
