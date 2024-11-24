classdef PPlane2d
    % PPLANE2D
    %
    % A class to plot phase portraits and trajectories for a system of
    % two-dimensional ordinary differential equations (ODEs) using variables u and v.
    %
    % The class accepts either symbolic expressions or function handles for du_dt and dv_dt.
    %
    % Examples:
    %   % Using symbolic expressions
    %   syms u v
    %   du_dt = -u + v;
    %   dv_dt = -2*v;
    %   pplane = PPlane2d(1, 2, 10, du_dt, dv_dt);
    %   pplane.plotAll();
    %
    %   % Using function handles
    %   du_dt = @(u, v) -u + v;
    %   dv_dt = @(u, v) -2*v;
    %   pplane = PPlane2d(1, 2, 10, du_dt, dv_dt);
    %   pplane.plotAll();

    properties
        u0          % Initial condition for u
        v0          % Initial condition for v
        Tend        % End time for ODE solver
        du_func     % Function handle for du/dt
        dv_func     % Function handle for dv/dt
        du_sym      % Symbolic expression for du/dt (if provided)
        dv_sym      % Symbolic expression for dv/dt (if provided)
    end

    properties (Hidden)
        uu          % Grid points for u
        vv          % Grid points for v
        duu         % Slope field values for du/dt
        dvv         % Slope field values for dv/dt
        tode        % Time vector from ODE solver
        yode        % Solution matrix from ODE solver
    end

    methods
        % CONSTRUCTOR
        function self = PPlane2d(u0, v0, Tend, du_dt, dv_dt)
            % PPLANE2D Constructor method
            %
            % Inputs:
            %   u0     - Initial condition for u
            %   v0     - Initial condition for v
            %   Tend   - End time for ODE solver
            %   du_dt  - Symbolic expression or function handle for du/dt
            %   dv_dt  - Symbolic expression or function handle for dv/dt

            if nargin ~= 5
                error('PPlane2d requires five input arguments: u0, v0, Tend, du_dt, dv_dt');
            end

            self.u0 = u0;
            self.v0 = v0;
            self.Tend = Tend;

            if isa(du_dt, 'function_handle') && isa(dv_dt, 'function_handle')
                % Inputs are function handles
                self.du_func = du_dt;
                self.dv_func = dv_dt;
            elseif isa(du_dt, 'sym') && isa(dv_dt, 'sym')
                % Inputs are symbolic expressions
                self.du_sym = du_dt;
                self.dv_sym = dv_dt;
                % Convert symbolic expressions to function handles
                try
                    self.du_func = matlabFunction(self.du_sym, 'Vars', {'u', 'v'});
                    self.dv_func = matlabFunction(self.dv_sym, 'Vars', {'u', 'v'});
                catch ME
                    error('Error converting symbolic expressions to function handles: %s', ME.message);
                end
            else
                error('du_dt and dv_dt must be both function handles or both symbolic expressions');
            end

            % Initialize grid and solve ODE
            self = self.initializeGrid();
            self = self.solveODE();
        end

        % INITIALIZE GRID
        function self = initializeGrid(self)
            % INITIALIZEGRID Initializes the grid for the phase portrait

            % Define grid range and density
            self.uu = -5:0.5:5;
            self.vv = -5:0.5:5;
            [self.uu, self.vv] = meshgrid(self.uu, self.vv);

            % Calculate slope fields on the grid
            try
                self.duu = self.du_func(self.uu, self.vv);
                self.dvv = self.dv_func(self.uu, self.vv);
            catch ME
                error('Error evaluating slope fields on the grid: %s', ME.message);
            end
        end

        % SOLVE ODE
        function self = solveODE(self)
            % SOLVEODE Solves the system of ODEs and stores the solution

            % Define the ODE system as a function handle
            ode_system = @(t, Y) [self.du_func(Y(1), Y(2)); self.dv_func(Y(1), Y(2))];

            % Initial conditions
            Y0 = [self.u0; self.v0];

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
            quiver(self.uu, self.vv, self.duu, self.dvv, 'r');
            hold on;
            xlabel('u');
            ylabel('v');
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
            plot(self.u0, self.v0, 'ko', 'MarkerFaceColor', 'g');  % Initial condition
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
                ylabel('u(t)');
                title('Solution Component u(t) over Time');
                grid on;
                
                subplot(2,1,2);
                plot(self.tode, self.yode(:,2), 'r+-', 'LineWidth', 1.5, 'MarkerSize', 5);
                xlabel('Time');
                ylabel('v(t)');
                title('Solution Component v(t) over Time');
                grid on;
                legend('v(t)');
            else
                % Overlay plots on the same figure with lines and markers
                figure;
                plot(self.tode, self.yode(:,1), 'ko-', 'LineWidth', 1.5, 'MarkerSize', 5);  % u(t) with circles
                hold on;
                plot(self.tode, self.yode(:,2), 'r+-', 'LineWidth', 1.5, 'MarkerSize', 5);  % v(t) with pluses
                xlabel('Time');
                ylabel('Solution Components');
                title('Solution Components over Time');
                legend('u(t)', 'v(t)');
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
            quiver(self.uu, self.vv, self.duu, self.dvv, 'r');
            hold on;
            xlabel('u');
            ylabel('v');
            title('Animated Phase Portrait with Evolving Trajectory');
            axis equal;
            grid on;
            
            % Initialize trajectory plot
            trajPlot = plot(NaN, NaN, 'b-', 'LineWidth', 2);
            % Plot initial condition
            initPlot = plot(self.u0, self.v0, 'ko', 'MarkerFaceColor', 'g');
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
            % ANIMATESOLUTIONOVERTIME Animates how the solutions u(t) and v(t) change over time
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
            
            % Determine fixed axes limits based on all data
            uLimits = [min(self.yode(:,1)), max(self.yode(:,1))];
            vLimits = [min(self.yode(:,2)), max(self.yode(:,2))];
            
            % Create figure for animation
            figure;
            subplot(2,1,1);
            xlabel('Time');
            ylabel('u(t)');
            title('Animated Evolution of u(t)');
            grid on;
            hold on;
            xlim([self.tode(1), self.tode(end)]);
            ylim(uLimits);
            plotU = plot(NaN, NaN, 'ko-', 'LineWidth', 1.5, 'MarkerSize', 5);
            
            subplot(2,1,2);
            xlabel('Time');
            ylabel('v(t)');
            title('Animated Evolution of v(t)');
            grid on;
            hold on;
            xlim([self.tode(1), self.tode(end)]);
            ylim(vLimits);
            plotV = plot(NaN, NaN, 'r+-', 'LineWidth', 1.5, 'MarkerSize', 5);
            
            % Animate u(t) and v(t)
            for k = 1:numSteps
                subplot(2,1,1);
                set(plotU, 'XData', self.tode(1:k), 'YData', self.yode(1:k,1));
                drawnow;
                pause(pauseTime);
                
                subplot(2,1,2);
                set(plotV, 'XData', self.tode(1:k), 'YData', self.yode(1:k,2));
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
            quiver(self.uu, self.vv, self.duu, self.dvv, 'r');
            hold on;
            xlabel('u');
            ylabel('v');
            title('Phase Portrait with Trajectory');
            axis equal;
            grid on;
            
            % Initialize trajectory plot
            trajPlot = plot(NaN, NaN, 'b-', 'LineWidth', 2);
            % Plot initial condition
            initPlot = plot(self.u0, self.v0, 'ko', 'MarkerFaceColor', 'g');

            legend('Slope Field', 'Trajectory', 'Initial Condition');
            
            % Right subplot: Solution components over time with animation
            subplot(1,2,2);
            xlabel('Time');
            ylabel('Solution Components');
            title('Evolution of Solution');
            grid on;
            hold on;
            
            % Determine fixed axes limits based on all data
            uLimits = [min(self.yode(:,1)), max(self.yode(:,1))];
            vLimits = [min(self.yode(:,2)), max(self.yode(:,2))];
            
            xlim([self.tode(1), self.tode(end)]);
            ylim([min([uLimits, vLimits]), max([uLimits, vLimits])]);  % Unified limits
            
            plotU = plot(NaN, NaN, 'ko-', 'LineWidth', 1.5, 'MarkerSize', 5);
            plotV = plot(NaN, NaN, 'r+-', 'LineWidth', 1.5, 'MarkerSize', 5);
            legend('u(t)', 'v(t)');
            
            % Animate both subplots
            for k = 1:numSteps
                % Update trajectory in phase portrait
                subplot(1,2,1);
                set(trajPlot, 'XData', self.yode(1:k,1), 'YData', self.yode(1:k,2));
                drawnow;
                
                % Update solution components over time
                subplot(1,2,2);
                set(plotU, 'XData', self.tode(1:k), 'YData', self.yode(1:k,1));
                set(plotV, 'XData', self.tode(1:k), 'YData', self.yode(1:k,2));
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
            quiver(self.uu, self.vv, self.duu, self.dvv, 'r');
            hold on;
            % Plot trajectory
            plot(self.yode(:,1), self.yode(:,2), 'b-', 'LineWidth', 2);
            % Plot initial condition
            plot(self.u0, self.v0, 'ko', 'MarkerFaceColor', 'g');
            xlabel('u');
            ylabel('v');
            title('Phase Portrait with Trajectory');
            axis equal;
            grid on;
            legend('Slope Field', 'Trajectory', 'Initial Condition');
            hold off;
            
            % Right subplot: Solution components over time
            subplot(1,2,2);
            % Plot solution components with lines and markers
            plot(self.tode, self.yode(:,1), 'ko-', 'LineWidth', 1.5, 'MarkerSize', 5);  % u(t) with circles
            hold on;
            plot(self.tode, self.yode(:,2), 'r+-', 'LineWidth', 1.5, 'MarkerSize', 5);  % v(t) with pluses
            xlabel('Time');
            ylabel('Solution Components');
            title('Solution Components over Time');
            legend('u(t)', 'v(t)');
            grid on;
            hold off;
        end
    end
end