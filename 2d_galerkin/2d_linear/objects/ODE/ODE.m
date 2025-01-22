classdef ODE

    properties
        g
        vInit
        s
        t0
        T
        timestep % current timestep
        constraints
        dt
        tGrid
        order
        solution
    end

    methods
        % CONSTRUCTOR
        function self = ODE(data,options)

            % store data
            self.g           = data.g;
            self.vInit       = data.vInit;
            self.s           = data.cofs.s;
            self.T           = data.time.T;
            self.dt          = data.time.dt;
            self.constraints = data.constraints;
            self.order = options.order;

            % compute additional data
            self.t0 = 0;
            self.tGrid = self.t0:self.dt:self.T;
            if self.tGrid(end) ~= self.T
                self.tGrid = [self.tGrid, self.T];
            end

            % initialize solution storage
            self.solution = zeros(1,length(self.tGrid));

        end

        function [t, u] = solve(obj)

            t = self.t0:self.dt:self.T;
            u = zeros(length(t), length(self.u0));
            u(1, :) = self.u0;
            for i = 2:length(t)
                u(i, :) = self.step(t(i), u(i-1, :));
            end

        end

        function V_out = solveIteration(self,V_in,SU)

				% remaining coefficients
				dt = self.dt;
				s  = self.s(SU,V_in);
				g  = self.g;
                V_prev = self.solution(self.timestep-1);

				% solve current iteration
				V_out = (dt * g + dt * s * SU + V_prev); 

				% apply resolvent
                vLower = self.constraints.vLower;
				vUpper = self.constraints.vUpper;
				V_out  = self.resolvent(V_out, dt, s, vLower, vUpper);

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
            indSkip = ceil(length(t) / 10);
			h = plot(t, v, '-', 'LineWidth', 4, 'MarkerSize', 10, 'MarkerIndices',[1:indSkip:length(t)]);

            % set x and y limits
            xMin = min(self.tGrid);
            xMax = max(self.tGrid);
            yMin = min(self.solution);
            yMax = max(self.solution);
            if yMin == yMax, yMax = yMin + 1; end
            xlim([xMin, xMax]);
			ylim([yMin - 0.1*yMin, yMax + 0.1*yMax]);

			% Label the axes and set a title
			xlabel('$t$','Interpreter','latex');
			ylabel('Temp.');
			title('$v(t)$','Interpreter','latex','FontSize',60);
			grid on;

		end

    end

    methods (Static)

        function v = resolvent(v_in,dt,c_B,vLower,vUpper)
        %RESOLVENT(V_IN,DT,C_B,VLOWER,VUPPER) computes the resolvent for the ODE
        %
        %			v' + c_b * (v - S(u)) + lambda = g
        %
        %	where lambda is a Lagrange multiplier enforcing some constraint.
        %
        % Author: Tyler Fara			Date: April 3, 2024
        %-----------------------------------------------------------------------------%
        % Inputs
        %	v_in	double or vector of doubles
        %	dt		double, represents time step size
        %	c_B		double, nonnegative
        %	vLower	double or NaN, lower constraint, i.e. v >= vLower
        %	vUpper  double or NaN, upper constraint, i.e. v <= vUpper
        %
        % Outputs
        %	v		double or vector of doubles
        %-----------------------------------------------------------------------------%
        % Notes
        %	(1) RESOLVENT is a general function that can be turned into a specific,
        %	anonymous function R(v_in) by choosing DT, C_B, VLOWER, VUPPER in advance.
        %	Supposing these variables have been assigned in advance, then the syntax is
        %
        %		R = @(v)(resolvent(v,dt,c_B,vLower,vUpper);
        %
        %	In this way, you can create custom resolvents enforcing lower, upper, or
        %	both kinds of constraints. And, these custom anonymous functions can be
        %	passed as arguments to other functions.
        %-----------------------------------------------------------------------------%

            denom = 1 + dt * c_B;

            % CASE 1: with upper constraint, no lower constraint
            if isnan(vLower) && ~isnan(vUpper)
                for i = 1:length(v_in)
                    if v_in(i) < vUpper * denom
                        v(i) = v_in(i) / denom;
                    else
                        v(i) = vUpper;
                    end
                end

            % CASE 2: no upper constraint, with lower constraint
            elseif ~isnan(vLower) && isnan(vUpper)
                for i = 1:length(v_in)
                    if v_in(i) > vLower * denom
                        v(i) = v_in(i) / denom;
                    else
                        v(i) = vLower;
                    end
                end

            % CASE 3: no upper constraint, no lower constraint
            elseif isnan(vLower) && isnan(vUpper)
                for i = 1:length(v_in)
                    v(i) = v_in(i) / denom;
                end

            % CASE 4: with upper constraint, with lower constraint
            elseif ~isnan(vLower) && ~isnan(vUpper)

                % check that constraints are appropriate
                if ~(vLower< vUpper)
                    error('vLower must be less than vUpper')
                end

                lowerBound = vLower * denom;
                upperBound = vUpper * denom;

                % apply resolvent
                for i = 1:length(v_in)
                    if v_in(i) > vLower * denom && v_in(i) < vUpper * denom
                        v(i) = v_in(i) / denom;
                    elseif v_in(i) <= vLower * denom
                        v(i) = vLower;
                    elseif v_in(i) >= vUpper * denom
                        v(i) = vUpper;
                    end
                end
            end
        end
    end
end