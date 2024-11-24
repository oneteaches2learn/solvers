% pplane_demo7.m
% Demonstration of the PPlane2d class with a specific ODE system:
% Incorporating r(u) and s(u) as:
%   r(u) = 0,    if u <= 0
%   r(u) = u,    if 0 < u < 1
%   r(u) = 1,    if u >= 1
% Similarly for s(u)

% Clear workspace and command window
clear;
clc;

% Define constants
k = 1;      % Diffusion constant
f = 0;      % Constant term in u equation
g = 1;      % Constant term in v equation

% Initial conditions
u0 = 1;     % Initial condition for u
v0 = 1;     % Initial condition for v

% Boundary Condition
u_L = 0;

% End time for ODE solver
Tend = 10;  % Adjust as needed

% Define r(u) and s(u) as function handles
r = @(u) (u <= 0).*0 + (u > 0 & u < 1).*u + (u >= 1).*1;
s = @(u) (u <= 0).*0 + (u > 0 & u < 1).*u + (u >= 1).*1;

% Define the ODEs as function handles using r(u) and s(u)
du_func = @(u, v) 3*f - (6*k).*(-u_L + 2*u - v) - (3*r(u)/12).*(u_L + 4*u + v) + (3*r(u)/2).*v;
dv_func = @(u, v) g - s(u).*v + (s(u)/4).*(u_L + 2*u + v);

% Create an instance of PPlane2d with the defined parameters
pplane = PPlane2d(u0, v0, Tend, du_func, dv_func);

% Plotting functions
% You can uncomment any of these as needed
% pplane.plotSlopeField();
% pplane.plotTrajectory();
% pplane.plotSolutionOverTime();
% pplane.plotSolutionOverTime(true);
pplane.plotAll();

% Animations
% pplane.animateTrajectory();
% pplane.animateSolutionOverTime();
% pplane.animateAll();