classdef PPlane2d
    % PPPLANE2D
    % 
    % A class to plot phase portraits and trajectories for a system of
    % two-dimensional ordinary differential equations (ODEs).
    %
    % Example:
    %   syms x y
    %   f = -x + y;
    %   g = -2*y;
    %   pplane = PPlane2d(1, 2, 10, f, g);
    %   pplane.solveODE();
    %   pplane.plotSlopeField();
    %   pplane.plotTrajectory();
    %   pplane.plotSolutionOverTime();

    properties
        x0          % Initial condition for x
        y0          % Initial condition for y
        Tend        % End time for ODE solver
        dx_sym      % Symbolic expression for dx/dt
        dy_sym      % Symbolic expression for dy/dt
        dx_func     % Function handle for dx/dt
        dy_func     % Function handle for dy/dt
        xx          % Grid points for x
        yy          % Grid points for y
        dxx         % Slope field values for dx/dt
        dyy         % Slope field values for dy/dt
        tode        % Time vector from ODE solver
        yode        % Solution matrix from ODE solver
    end
    
    properties (Hidden)
        %{
        dx_func     % Function handle for dx/dt
        dy_func     % Function handle for dy/dt
        xx          % Grid points for x
        yy          % Grid points for y
        dxx         % Slope field values for dx/dt
        dyy         % Slope field values for dy/dt
        tode        % Time vector from ODE solver
        yode        % Solution matrix from ODE solver
        %}
    end
    
    methods
        % CONSTRUCTOR
        function self = PPlane2d(x0, y0, Tend, dx_sym, dy_sym)
            % PPPLANE2D Constructor method
            %
            % Inputs:
            %   x0     - Initial condition for x
            %   y0     - Initial condition for y
            %   Tend   - End time for ODE solver
            %   dx_sym - Symbolic expression for dx/dt
            %   dy_sym - Symbolic expression for dy/dt
            
            if nargin ~= 5
                error('PPlane2d requires five input arguments: x0, y0, Tend, dx_sym, dy_sym');
            end
            
            self.x0 = x0;
            self.y0 = y0;
            self.Tend = Tend;
            self.dx_sym = dx_sym;
            self.dy_sym = dy_sym;
            
            % Convert symbolic expressions to function handles
            try
                self.dx_func = matlabFunction(self.dx_sym, 'Vars', {'x', 'y'});
                self.dy_func = matlabFunction(self.dy_sym, 'Vars', {'x', 'y'});
            catch ME
                error('Error converting symbolic expressions to function handles: %s', ME.message);
            end
            
            % Initialize grid
            self = self.initializeGrid();
            self = self.solveODE();
        end
        
        % INITIALIZE GRID
        function self = initializeGrid(self)
            % Initialize the grid for the phase portrait
            self.xx = -5:0.5:5;
            self.yy = -5:0.5:5;
            [self.xx, self.yy] = meshgrid(self.xx, self.yy);
            
            % Calculate slope fields on the grid
            try
                self.dxx = self.dx_func(self.xx, self.yy);
                self.dyy = self.dy_func(self.xx, self.yy);
            catch ME
                error('Error evaluating slope fields on the grid: %s', ME.message);
            end
        end
        
        % SOLVE ODE
        function self = solveODE(self)
            % SOLVEODE Solves the system of ODEs and stores the solution
            
            % Define the ODE system as a function handle
            ode_system = @(t, Y) [self.dx_func(Y(1), Y(2)); self.dy_func(Y(1), Y(2))];
            
            % Initial conditions
            Y0 = [self.x0; self.y0];
            
            % Solve the ODE using ode45
            try
                [self.tode, self.yode] = ode45(ode_system, [0, self.Tend], Y0);
            catch ME
                error('Error solving ODEs: %s', ME.message);
            end
        end
        
        % PLOT SLOPE FIELD
        function plotSlopeField(self)
            % PLOTSLOPEFIELD Plots the slope field using quiver
            figure;
            quiver(self.xx, self.yy, self.dxx, self.dyy, 'r');
            hold on;
            xlabel('x');
            ylabel('y');
            title('Phase Portrait');
            axis equal;
            grid on;
        end
        
        % PLOT TRAJECTORY
        function plotTrajectory(self)
            % PLOTTRAJECTORY Plots the trajectory of the ODE solution on the phase portrait
            if isempty(self.yode)
                error('ODE not solved yet. Call solveODE() before plotting trajectory.');
            end
            
            plot(self.yode(:,1), self.yode(:,2), 'b-', 'LineWidth', 2);
            plot(self.x0, self.y0, 'ko', 'MarkerFaceColor', 'g');  % Initial condition
            legend('Slope Field', 'Trajectory', 'Initial Condition');
            hold off;
        end
        
        % PLOT SOLUTION OVER TIME
        function plotSolutionOverTime(self)
            % PLOTSOLUTIONOVERTIME Plots the solution components over time
            if isempty(self.yode)
                error('ODE not solved yet. Call solveODE() before plotting solutions.');
            end
            
            figure;
            subplot(2,1,1);
            plot(self.tode, self.yode(:,1), 'k-', 'LineWidth', 1.5);
            xlabel('Time');
            ylabel('x(t)');
            title('Solution Components over Time');
            grid on;
            
            subplot(2,1,2);
            plot(self.tode, self.yode(:,2), 'r+', 'LineWidth', 1.5);
            xlabel('Time');
            ylabel('y(t)');
            title('Solution Components over Time');
            grid on;
            legend('y(t)');
        end
        
        % PLOT ALL
        function plotAll(self)
            % PLOTALL Executes all plotting methods
            self.plotSlopeField();
            if isempty(self.yode)
                self.solveODE();
            end
            self.plotTrajectory();
            self.plotSolutionOverTime();
        end
    end
end