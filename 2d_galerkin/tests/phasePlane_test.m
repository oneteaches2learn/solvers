% Define symbolic variables
syms x y

% Define the ODEs symbolically
f = -x + y;    % Example: x' = -x + y
g = -2*y;      % Example: y' = -2y

% Initial conditions and end time
x0 = 1;
y0 = 2;
Tend = 10;

% Call the modified my_pplane function
my_pplane(x0, y0, Tend, f, g);