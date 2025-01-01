classdef ODE

    properties
        g
        vInit
        s
        t0
        T
        dt
        tGrid
        order
        solution
    end

    methods
        % CONSTRUCTOR
        function self = ODE(data,options)

            % store data
            self.g = data.g;
            self.vInit = data.vInit;
            self.s = data.cofs.s;
            self.T = data.time.T;
            self.dt = data.time.dt;
            self.order = options.order;

            % compute additional data
            self.t0 = 0;
            self.tGrid = self.t0:self.dt:self.T;
            if self.tGrid(end) ~= self.T
                self.tGrid = [self.tGrid, self.T];
            end

        end

        function [t, u] = solve(obj)

            t = self.t0:self.dt:self.T;
            u = zeros(length(t), length(self.u0));
            u(1, :) = self.u0;
            for i = 2:length(t)
                u(i, :) = self.step(t(i), u(i-1, :));
            end

        end

        function u = step(obj, t, u)

            switch obj.order
                case 1
                    u = u + obj.dt * obj.g(t, u);
                case 2
                    u = u + obj.dt * obj.g(t, u) + 0.5 * obj.dt^2 * obj.g(t, u);
                otherwise
                    error('Invalid order');
            end

        end

        % SETTERS
        function self = set.tGrid(self,tGrid)
        % set.tGrid allows to pass a time grid
        %    If the ODE is coupled to a PDE, a time grid is passed from the PDE
        %    to ensure no discrepancies between the time grids used in the PDE
        %    and ODE.
        
            self.tGrid = tGrid;
        end


        % PLOTTING FUNCTIONS
		function plot(self,timestep)
		% PLOTODE Plots the solution to the ODE

            % Set default timestep to the end of the simulation
            if nargin == 1
                timestep = length(self.tGrid);
            end

            % store data
            t = self.tGrid(1:timestep);
            v = self.solution(1:timestep);

			% Plot self.V on the y-axis against time on the x-axis
			plot(t, v, 'LineWidth', 2);

            % set x and y limits
            xMin = min(self.tGrid);
            xMax = max(self.tGrid);
            yMin = min(self.solution);
            yMax = max(self.solution);
            if xMin == xMax, xMax = xMin + 1; end
            if yMin == yMax, yMax = yMin + 1; end
            xlim([xMin, xMax]);
            ylim([yMin, yMax]);

			% Label the axes and set a title
			xlabel('Time');
			ylabel('V');
			title('v vs. Time');
			grid on;

		end

    end
end