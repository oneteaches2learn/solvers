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
    %   pplane = pplane.solveODE();
    %   pplane.plotSlopeField();
    %   pplane.plotTrajectory();
    %   pplane.plotSolutionOverTime();
    %
    %   % Plot everything using subplots
    %   pplane.plotAll();
    %
    %   % Animate the trajectory
    %   pplane.animateTrajectory();
    %
    %   % Animate the solution over time
    %   pplane.animateSolutionOverTime();
    %
    %   % Animate all plots together
    %   pplane.animateAll();

    properties
        x0          % Initial condition for x
        y0          % Initial condition for y
        Tend        % End time for ODE solver
        dx_sym      % Symbolic expression for dx/dt
        dy_sym      % Symbolic expression for dy/dt
    end
    
    properties (Hidden)
        dx_func     % Function handle for dx/dt
        dy_func     % Function handle for dy/dt
        xx          % Grid points for x
        yy          % Grid points for y
        dxx         % Slope field values for dx/dt
        dyy         % Slope field values for dy/dt
        tode        % Time vector from ODE solver
        yode        % Solution matrix from ODE solver
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
            
            % Initialize grid and solve ODE
            self = self.initializeGrid();
            self = self.solveODE();
        end
        
        % INITIALIZE GRID
        function self = initializeGrid(self)
            % INITIALIZEGRID Initializes the grid for the phase portrait
            %
            % Outputs:
            %   self - Updated object with grid and slope fields
            
            % Define grid range and density
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
            %
            % Outputs:
            %   self - Updated object with ODE solution
            
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
        function plotSolutionOverTime(self, useSubplots)
            % PLOTSOLUTIONOVERTIME Plots the solution components over time
            %
            % Usage:
            %   plotSolutionOverTime()             % Overlay plots (default)
            %   plotSolutionOverTime(true)         % Use subplots
            
            if nargin < 2
                useSubplots = false;  % Default behavior: overlay
            end
            
            if isempty(self.yode)
                error('ODE not solved yet. Call solveODE() before plotting solutions.');
            end
            
            if useSubplots
                % Display plots separately using subplots
                figure;
                
                subplot(2,1,1);
                plot(self.tode, self.yode(:,1), 'ko-', 'LineWidth', 1.5, 'MarkerSize', 5);
                xlabel('Time');
                ylabel('x(t)');
                title('Solution Component x(t) over Time');
                grid on;
                
                subplot(2,1,2);
                plot(self.tode, self.yode(:,2), 'r+-', 'LineWidth', 1.5, 'MarkerSize', 5);
                xlabel('Time');
                ylabel('y(t)');
                title('Solution Component y(t) over Time');
                grid on;
                legend('y(t)');
            else
                % Overlay plots on the same figure with lines and markers
                figure;
                plot(self.tode, self.yode(:,1), 'ko-', 'LineWidth', 1.5, 'MarkerSize', 5);  % x(t) with circles
                hold on;
                plot(self.tode, self.yode(:,2), 'r+-', 'LineWidth', 1.5, 'MarkerSize', 5);  % y(t) with pluses
                xlabel('Time');
                ylabel('Solution Components');
                title('Solution Components over Time');
                legend('x(t)', 'y(t)');
                grid on;
                hold off;
            end
        end
        
        % ANIMATE TRAJECTORY
        function animateTrajectory(self)
            % ANIMATETRAJECTORY Animates the trajectory evolving through the slope field over time
            %
            % Usage:
            %   animateTrajectory()
            
            if isempty(self.yode)
                error('ODE not solved yet. Call solveODE() before animating trajectory.');
            end
            
            % Total animation time in seconds
            totalTime = 1;
            numSteps = length(self.tode);
            pauseTime = totalTime / numSteps;
            
            % Create figure for animation
            figure;
            quiver(self.xx, self.yy, self.dxx, self.dyy, 'r');
            hold on;
            xlabel('x');
            ylabel('y');
            title('Animated Phase Portrait with Evolving Trajectory');
            axis equal;
            grid on;
            
            % Initialize trajectory plot
            trajPlot = plot(NaN, NaN, 'b-', 'LineWidth', 2);
            % Plot initial condition
            initPlot = plot(self.x0, self.y0, 'ko', 'MarkerFaceColor', 'g');
            legend('Slope Field', 'Trajectory', 'Initial Condition');
            
            % Animate trajectory
            for k = 1:numSteps
                set(trajPlot, 'XData', self.yode(1:k,1), 'YData', self.yode(1:k,2));
                drawnow;
                pause(pauseTime);
            end
            hold off;
        end
        
        % ANIMATE SOLUTION OVER TIME
        function animateSolutionOverTime(self)
            % ANIMATESOLUTIONOVERTIME Animates how the solutions x(t) and y(t) change over time
            %
            % Usage:
            %   animateSolutionOverTime()
            
            if isempty(self.yode)
                error('ODE not solved yet. Call solveODE() before animating solutions.');
            end
            
            % Total animation time in seconds
            totalTime = 1;
            numSteps = length(self.tode);
            pauseTime = totalTime / numSteps;
            
            % Determine fixed axes limits based on the first time step
            xLimits = [min(self.yode(:,1)), max(self.yode(:,1))];
            yLimits = [min(self.yode(:,2)), max(self.yode(:,2))];
            
            % Create figure for animation
            figure;
            subplot(2,1,1);
            xlabel('Time');
            ylabel('x(t)');
            title('Animated Evolution of x(t)');
            grid on;
            hold on;
            xlim([self.tode(1), self.tode(end)]);
            ylim(xLimits);
            plotX = plot(NaN, NaN, 'ko-', 'LineWidth', 1.5, 'MarkerSize', 5);
            
            subplot(2,1,2);
            xlabel('Time');
            ylabel('y(t)');
            title('Animated Evolution of y(t)');
            grid on;
            hold on;
            xlim([self.tode(1), self.tode(end)]);
            ylim(yLimits);
            plotY = plot(NaN, NaN, 'r+-', 'LineWidth', 1.5, 'MarkerSize', 5);
            
            % Animate x(t)
            for k = 1:numSteps
                subplot(2,1,1);
                set(plotX, 'XData', self.tode(1:k), 'YData', self.yode(1:k,1));
                drawnow;
                pause(pauseTime);
                
                subplot(2,1,2);
                set(plotY, 'XData', self.tode(1:k), 'YData', self.yode(1:k,2));
                drawnow;
                pause(pauseTime);
            end
            hold off;
        end
        
        % ANIMATE ALL
        function animateAll(self)
            % ANIMATEALL Animates both the trajectory and the solutions over time
            %
            % Usage:
            %   animateAll()
            
            if isempty(self.yode)
                error('ODE not solved yet. Call solveODE() before animating.');
            end
            
            % Total animation time in seconds
            totalTime = 1;
            numSteps = length(self.tode);
            pauseTime = totalTime / numSteps;
            
            % Create figure with subplots
            figure('Position', [100, 100, 1200, 600]);  % [left, bottom, width, height]
            
            % Left subplot: Phase portrait with animation
            subplot(1,2,1);
            quiver(self.xx, self.yy, self.dxx, self.dyy, 'r');
            hold on;
            xlabel('x');
            ylabel('y');
            title('Phase Portrait with Trajectory');
            axis equal;
            grid on;
            
            % Initialize trajectory plot
            trajPlot = plot(NaN, NaN, 'b-', 'LineWidth', 2);
            % Plot initial condition
            initPlot = plot(self.x0, self.y0, 'ko', 'MarkerFaceColor', 'g');
            legend('Slope Field', 'Trajectory', 'Initial Condition');
            
            % Right subplot: Solution components over time with animation
            subplot(1,2,2);
            xlabel('Time');
            ylabel('Solution Components');
            title('Evolution of Solution');
            grid on;
            hold on;
            
            % Determine fixed axes limits based on the first time step
            xLimits = [min(self.yode(:,1)), max(self.yode(:,1))];
            yLimits = [min(self.yode(:,2)), max(self.yode(:,2))];
            
            xlim([self.tode(1), self.tode(end)]);
            ylim([min([xLimits, yLimits]), max([xLimits, yLimits])]);  % Unified limits
            
            plotX = plot(NaN, NaN, 'ko-', 'LineWidth', 1.5, 'MarkerSize', 5);
            plotY = plot(NaN, NaN, 'r+-', 'LineWidth', 1.5, 'MarkerSize', 5);
            legend('x(t)', 'y(t)');
            
            % Animate both subplots
            for k = 1:numSteps
                % Update trajectory in phase portrait
                subplot(1,2,1);
                set(trajPlot, 'XData', self.yode(1:k,1), 'YData', self.yode(1:k,2));
                drawnow;
                
                % Update solution components over time
                subplot(1,2,2);
                set(plotX, 'XData', self.tode(1:k), 'YData', self.yode(1:k,1));
                set(plotY, 'XData', self.tode(1:k), 'YData', self.yode(1:k,2));
                drawnow;
                
                pause(pauseTime);
            end
            hold off;
        end
        
        % PLOT ALL
        function plotAll(self)
            % PLOTALL Executes plotting trajectory and solution over time using subplots
            %
            % Displays:
            %   - Left subplot: Phase portrait with trajectory
            %   - Right subplot: Solution components over time
                
            if isempty(self.yode)
                self = self.solveODE();
            end
            
            % Create a new figure with specified location
            figure('Position', [100, 100, 1200, 600]);  % [left, bottom, width, height]
            
            % Left subplot: Phase portrait with trajectory
            subplot(1,2,1);
            % Plot slope field
            quiver(self.xx, self.yy, self.dxx, self.dyy, 'r');
            hold on;
            % Plot trajectory
            plot(self.yode(:,1), self.yode(:,2), 'b-', 'LineWidth', 2);
            % Plot initial condition
            plot(self.x0, self.y0, 'ko', 'MarkerFaceColor', 'g');
            xlabel('x');
            ylabel('y');
            title('Phase Portrait with Trajectory');
            axis equal;
            grid on;
            legend('Slope Field', 'Trajectory', 'Initial Condition');
            hold off;
            
            % Right subplot: Solution components over time
            subplot(1,2,2);
            % Plot solution components with lines and markers
            plot(self.tode, self.yode(:,1), 'ko-', 'LineWidth', 1.5, 'MarkerSize', 5);  % x(t) with circles
            hold on;
            plot(self.tode, self.yode(:,2), 'r+-', 'LineWidth', 1.5, 'MarkerSize', 5);  % y(t) with pluses
            xlabel('Time');
            ylabel('Solution Components');
            title('Solution Components over Time');
            legend('x(t)', 'y(t)');
            grid on;
            hold off;
        end
    end
end