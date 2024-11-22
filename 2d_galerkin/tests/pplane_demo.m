% Define symbolic variables
syms x y

% Define the ODEs symbolically
f = -x + y;    % Example: dx/dt = -x + y
g = -2*y;      % Example: dy/dt = -2y

% Initial conditions and end time
x0 = 1;
y0 = 2;
Tend = 10;

% Create an instance of PPlane2d
pplane = PPlane2d(x0, y0, Tend, f, g);

% Solve the ODE
pplane.solveODE();
fprintf('hi\n');

% Plot the slope field
%pplane.plotSlopeField();
%pplane.plotTrajectory();
pplane.plotSolutionOverTime(1);

% Alternatively, plot everything at once
pplane.plotAll();