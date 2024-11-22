function my_pplane(x0, y0, Tend, dx_sym, dy_sym)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % M. Peszynska for MTH 4/581
    % Modified to accept symbolic ODEs
    % plot phase portrait for x'=f(x,y) and y'=g(x,y) and a trajectory with some I.C.
    % example:
    %   syms f g
    %   f = -x;
    %   g = -2*y;
    %   my_pplane(1,2,5, f, g)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Validate input arguments
    if nargin ~= 5
        error('my_pplane requires five input arguments: x0, y0, Tend, dx_sym, dy_sym');
    end

    % Create a grid for the phase portrait
    [xx, yy] = meshgrid(-5:0.5:5, -5:0.5:5);  % Increased range for better visualization

    % Convert symbolic expressions to function handles
    % Assume that dx_sym and dy_sym are symbolic expressions involving x and y
    % Convert them to MATLAB function handles using matlabFunction
    try
        dx_func = matlabFunction(dx_sym, 'Vars', {'x', 'y'});
        dy_func = matlabFunction(dy_sym, 'Vars', {'x', 'y'});
    catch ME
        error('Error converting symbolic expressions to function handles: %s', ME.message);
    end

    % Calculate slope fields on the grid
    dxx = dx_func(xx, yy);
    dyy = dy_func(xx, yy);

    % Plot the slope fields using quiver
    figure;
    quiver(xx, yy, dxx, dyy, 'r');
    hold on;
    xlabel('x');
    ylabel('y');
    title('Phase Portrait and Trajectory');
    axis equal;
    grid on;

    % Pause to allow user to inspect the slope field
    fprintf('Slope field plotted. Press any key to continue with trajectory...\n');
    pause;

    % Define the system of ODEs as a function handle for ode45
    ode_system = @(t, Y) [dx_func(Y(1), Y(2)); dy_func(Y(1), Y(2))];

    % Solve the ODE using ode45
    try
        [tode, yode] = ode45(ode_system, [0, Tend], [x0; y0]);
    catch ME
        error('Error solving ODEs: %s', ME.message);
    end

    % Plot the trajectory on the phase portrait
    plot(yode(:,1), yode(:,2), 'b-', 'LineWidth', 2);
    plot(x0, y0, 'ko', 'MarkerFaceColor', 'g');  % Initial condition
    legend('Slope Field', 'Trajectory', 'Initial Condition');
    pause;

    % Plot the solutions over time
    figure;
    subplot(2,1,1);
    plot(tode, yode(:,1), 'k-', 'LineWidth', 1.5);
    xlabel('Time');
    ylabel('x(t)');
    title('Solution Components over Time');
    grid on;

    subplot(2,1,2);
    plot(tode, yode(:,2), 'r+', 'LineWidth', 1.5);
    xlabel('Time');
    ylabel('y(t)');
    grid on;
    legend('y(t)');
end