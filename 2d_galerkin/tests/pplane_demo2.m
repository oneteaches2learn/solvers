% pplane_demo2.m
% Demonstration of the PPlane2d class with a specific ODE system:
% dμ/dt + r(μ - v) + b(Lμ - v) = f
% dv/dt + s(v - μ) = g
%
% Where:
%   r > 0 and s > 0 are positive constants
%   f and g are constants
%
% The system is rewritten in standard form as:
%   dμ/dt = -r(μ - v) - b(Lμ - v) + f
%   dv/dt = -s(v - μ) + g

% Clear workspace and command window, define symbolic variables
clear; 
syms u v

% USER INPUTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define constants
r = 1;      % Positive constant for μ equation
s = 1;    % Positive constant for v equation
b = 1;
L = u;
f = 0;      % Constant term in μ equation
g = 0;      % Constant term in v equation

% Initial conditions
u0 = -2;    % Initial condition for μ
v0 = 2;    % Initial condition for v

% End time for ODE solver
Tend = 2;  % You can adjust this as needed


% BLACK BOX %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define the ODEs in standard form
du_dt = -r*(u - v) - b*(L*u - v) + f;
dv_dt = -s*(v - u) + g;

% Create an instance of PPlane2d with the defined parameters
pplane = PPlane2d(u0, v0, Tend, du_dt, dv_dt);


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