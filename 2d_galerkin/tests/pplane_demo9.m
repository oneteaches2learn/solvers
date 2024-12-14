% pplane_demo1.m
% Demonstration of the PPlane2d class with a specific ODE system:
% dμ/dt + |partial Omega|/|Omega| (alpha μ + beta v + gamma) + r(μ - v) = f
% dv/dt + s(v - μ) = g
%
% Where:
%   r > 0 and s > 0 are positive constants
%   f and g are constants
%
% The system is rewritten in standard form as:
%   dμ/dt = -r(μ - v) + f
%   dv/dt = -s(v - μ) + g

% Clear workspace and command window, define symbolic variables
clear; 
syms u v

% USER INPUTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define constants
r = 1;      % Positive constant for μ equation
s = 1;     % Positive constant for v equation
A = 1;      % A = |partial Omega|/|Omega|
alpha = 0;
beta  = 10;
gamma = 0;
f = 0;      % Constant term in μ equation
g = 1;      % Constant term in v equation

% Initial conditions
u0 = 2;    % Initial condition for μ
v0 = 2;    % Initial condition for v

% End time for ODE solver
Tend = 10;  % You can adjust this as needed


% BLACK BOX %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define the ODEs in standard form
du_dt = f - r*(u - v) - A*(alpha*u + beta*v + gamma);
dv_dt = g - s*(v - u);

% Create an instance of PPlane2d with the defined parameters
pplane = PPlane2d(u0, v0, Tend, du_dt, dv_dt);

% Compute Steady State
mu_s = (f - (A * beta * g / s + A * gamma - r * g / s)) / (A * (alpha + beta))
v_s = (g + s * mu_s) / s

% POST-PROCESSING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plotting functions
%pplane.plotSlopeField();
%pplane.plotTrajectory();
%pplane.plotSolutionOverTime();
%pplane.plotSolutionOverTime(true);
pplane.plotAll();

% Animations
%pplane.animateTrajectory();
%pplane.animateSolutionOverTime();
%pplane.animateAll();