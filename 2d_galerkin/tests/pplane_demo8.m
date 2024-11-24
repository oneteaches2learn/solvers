% pplane_demo8.m
% Demonstration of the PPlane2d class with a specific ODE system:
% Incorporating r(v) and s(v) as:
%   r(v) = 0,    if v <= 0
%   r(v) = v,    if 0 < v < 1
%   r(v) = 1,    if v >= 1
% Similarly for s(v)

% Clear workspace and command window
clear;
clc;

% Define constants
k = 1;      % Diffusion constant
f = 0;      % Constant term in u equation
g = 0;      % Constant term in v equation

% Initial conditions
u0 = 0;     % Initial condition for u
v0 = 0;     % Initial condition for v

% Boundary Condition
u_L = 0;

% End time for ODE solver
Tend = 10;  % Adjust as needed

% Define r(v) and s(v) as function handles
r = @(v) (v <= 0).*0 + (v > 0 & v < 1).*v + (v >= 1).*1;
s = @(v) (v <= 0).*0 + (v > 0 & v < 1).*v + (v >= 1).*1;

% Define the ODEs as function handles using r(v) and s(v)
du_func = @(u, v) 3*f - (6*k).*(-u_L + 2*u - v) - (3*r(v)/12).*(u_L + 4*u + v) + (3*r(v)/2).*v;
dv_func = @(u, v) g - s(v).*v + (s(v)/4).*(u_L + 2*u + v);

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