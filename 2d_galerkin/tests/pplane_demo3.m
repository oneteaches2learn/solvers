% pplane_demo3.m
% Demonstration of the PPlane2d class with a specific ODE system:
%   du/dt = 3*f - 6*( 2*k*u - k*u_L - k*v)
%   dv/dt = g - s*v + (s/4)*(u_L + 2*u + v)
%
% Where:
%   r > 0 and s > 0 are positive constants
%   f and g are constants
%
% This system represents the system 
%   u_t - nabla (k nabla u) = f %   v_t + s(u - S(u)) = g
% 
% with boundary conditions
%   u(0) = u0
%   u(L) = v
%
% on a 1D domain with 3 nodes.

clear; syms u v
% USER INPUTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define constants
r = 1;      % Positive constant for μ equation
s = 1;      % Positive constant for v equation
f = 0;      % Constant term in μ equation
g = 1;      % Constant term in v equation

% Initial conditions
u0 = 1;    % Initial condition for μ
v0 = 1;    % Initial condition for v

% Boundary Condition
u_L = -1;

% End time for ODE solver
Tend = 10;  % You can adjust this as needed


% BLACK BOX %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define the ODEs in standard form
du_dt = 3*f - 6*(2*r*u - r*u_L - r*v);
dv_dt = g - s*v + s/4 * (u_L + 2*u + v);

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