% pplane_demo.m
% Demonstration of the PPlane2d class with a specific ODE system:
% dμ/dt + r(μ - v) = f
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
syms mu v

% USER INPUTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define constants
r = 1;      % Positive constant for μ equation
s = 1.5;    % Positive constant for v equation
f = 1;      % Constant term in μ equation
g = 1;      % Constant term in v equation

% Initial conditions
mu0 = 1;    % Initial condition for μ
v0 = 1;     % Initial condition for v

% End time for ODE solver
Tend = 10;  % You can adjust this as needed


% BLACK BOX %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define the ODEs in standard form
dmu_dt = -r*(mu - v) + f;
dv_dt = -s*(v - mu) + g;

% Create an instance of PPlane2d with the defined parameters
pplane = PPlane2d(mu0, v0, Tend, dmu_dt, dv_dt);

%{
% Solve the ODE
pplane.solveODE();


% POST-PROCESSING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plotting functions
%pplane.plotSlopeField();
%pplane.plotTrajectory();
%pplane.plotSolutionOverTime();
%pplane.plotSolutionOverTime(true);
%pplane.plotAll();

% Animations
%pplane.animateTrajectory();
%pplane.animateSolutionOverTime();
pplane.animateAll();
%}