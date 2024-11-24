% pplane_demo4.m
% Demonstration of the PPlane2d class with a specific ODE system:
%   du/dt = 3*f - (6*k)*(-u_L + 2*u - v) - (3*r/12)*(u_L + 4*u + u_R) + (3*r/2)*v
%   dv/dt = g - s*v + (s/4)*(u_L + 2*u + u_R)
% 
%
% Where:
%   r > 0 and s > 0 are positive constants
%   f and g are constants
%
% This system represents the system 
%   u_t - nabla (k nabla u) + r(u - v) = f
%   v_t + s(u - S(u)) = g
% 
% with boundary conditions
%   u(0) = u_L
%   u(1) = u_R
%
% on a 1D domain with 3 nodes.

clear; syms u v
% USER INPUTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define constants
r = 1;      % Positive constant for μ equation
s = 1;    % Positive constant for v equation
f = 0;      % Constant term in μ equation
g = 1;      % Constant term in v equation

% Initial conditions
u0 = 2;    % Initial condition for μ
v0 = 2;    % Initial condition for v

% Boundary Condition
u_L = 0;
u_R = v;

% End time for ODE solver
Tend = 10;  % You can adjust this as needed


% BLACK BOX %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define the ODEs in standard form
du_dt = 3*f - (6*r)*(-u_L + 2*u - v) - (3*r/12)*(u_L + 4*u + u_R) + (3*r/2)*v;
dv_dt = g - s*v + (s/4)*(u_L + 2*u + u_R);

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