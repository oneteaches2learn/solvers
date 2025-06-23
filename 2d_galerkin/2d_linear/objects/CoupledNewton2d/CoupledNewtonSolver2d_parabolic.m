classdef CoupledNewtonSolver2d_parabolic < NewtonGalerkinSolver2d_parabolic

	properties
		PDE
		ODE
		V_prev % value of V at previous time step
		V	   % value of V at current time and current iteration
		s_history % <~~ temporary?
	end

	methods
		function self = CoupledNewtonSolver2d_parabolic(dom,auxfun,ODE)
			
			% call superclass constructor
			self@NewtonGalerkinSolver2d_parabolic(dom,auxfun);

			% store additional data
			self.ODE = ODE;
			self.ODE.tGrid = self.domain.time.tGrid;

			% solve
			% ... handled in individual subclasses

		end

        function self = solveTimestep(self,S,b,FreeNodes)

            % store dirichlet nodes and free nodes
			FreeNodes = self.domain.boundary.freeNodes;
			dirichlet = unique(self.domain.boundary.D_nodes);

            % construct U_tilde
            U_tilde = self.U;
            U_tilde(dirichlet) = 0;

            % Newton-Galerkin loop
            for iter = 1:100

				% solve PDE iteration
                [self, U_tilde, W] = self.solveIteration(U_tilde);

				% solve ODE iteration
				SU = self.domain.nodalQuadrature(self.U) / self.domain.domainArea;
				self.V = self.ODE.solveIteration(self.V,SU,self.t);

                % check convegence
                if norm(W) < 10^(-10)
                    %fprintf(' timestep: %d, converged in %d iterations\n',self.timestep,iter)
                    self.iterHistory(self.timestep) = iter;
                    break
                end
            end

            % store resolved solution
            self.solution(:,self.timestep) = self.U;
			self.ODE.solution(self.timestep) = self.V;

        end

		function self = initializeProblem(self)

			% initialize PDE by calling superclass method
			self = self.initializeProblem@NewtonGalerkinSolver2d_parabolic();

			% initialize ODE 
			self.V_prev = self.ODE.vInit;
			self.V = self.V_prev;
			self.ODE.solution(1) = self.V;
			self.ODE.timestep = self.timestep;

		end

		function self = initializeTimestep(self)

			% initialize PDE timestep by calling superclass method
			self = self.initializeTimestep@NewtonGalerkinSolver2d_parabolic();

			% initialize ODE timestep
			self.V_prev = self.V;
			self.ODE.timestep = self.timestep;

		end

		% PLOTTING FUNCTIONS
		function plot(self,timestep)

			% default to last timestep
			if nargin == 1
				timestep = self.domain.time.N_t;
			end

			% call superclass method
			plot@NewtonGalerkinSolver2d_parabolic(self,timestep);

			% Adjust figure
			view(2);

			% remove tick labels
            set(gca,'YTickLabel',[]);
            set(gca,'XTickLabel',[]);
			
			% adjust title
			str = sprintf('$u(x,t), t = %.0f$',self.domain.time.tGrid(timestep));
			title(str,'Interpreter','latex','FontSize',80);

		end

		function plotPDE(self,timestep)

			% default to last timestep
			if nargin == 1
				timestep = self.domain.time.N_t;
			end

			% call plot method
			self.plot(timestep);

		end

		function plotODE(self,timestep)

			% Default to the last timestep
			if nargin == 1
				timestep = self.domain.time.N_t;
			end

			self.ODE.plot(timestep);

		end

		function plotTemperatures(self,timestep)

			% Default to the last timestep
			if nargin == 1
				timestep = self.domain.time.N_t;
			end

			% combine plots
			hold on
			self.plotODE(timestep);
			h2 = self.plotAverageSolutionValue(timestep);
			hold off

			% change second plot's linestyle
			h2.LineStyle = '--';

			% set title, legend, and labels
			title('$v$ and $\langle u \rangle_\Omega$','Interpreter','latex','FontSize',50);
			legend('v','$\langle u \rangle_\Omega$','Interpreter','latex','Location','east');

			% set ylim
			yMin = min(min(self.ODE.solution),min(self.getAverageSolutionValue));
			yMax = max(max(self.ODE.solution),max(self.getAverageSolutionValue));
			ylim([yMin - 0.1 * (yMax - yMin), yMax + 0.1 * (yMax - yMin)]);

			f = gcf;
			f.Position = [100, 100, 400, 400];
			%Legend = f.Children(1);
			%Legend.Location = 'best';

		end

		function plotBoth(self,timestep)
		% ANIMATE Sequentially plots the PDE + ODE solutions from the first
		% to the last timestep. The total animation duration is one second.
				
			mode.plot1 = 'PDE';
			mode.plot2 = 'Temperatures';

			% set first plot function
			if strcmp(mode.plot1,'PDE')
				plot1 = @(t) self.plotPDE(t);
			elseif strcmp(mode.plot1,'Frostbite')
				plot1 = @(t) self.plotFrostbite(t);
			elseif strcmp(mode.plot1,'BloodSupply')
				plot1 = @(t) self.plotBloodSupply(t);
			end

			% set second plot function
			if strcmp(mode.plot2,'ODE')
				plot2 = @(t) self.plotODE(t);
			elseif strcmp(mode.plot2,'Temperatures')
				plot2 = @(t) self.plotTemperatures(t);
			elseif strcmp(mode.plot2,'Average')
				plot2 = @(t) self.plotAverageSolutionValue(t);
			end

			% Number of timesteps
			nT = self.domain.time.N_t;
			
			% Total animation time (in seconds)
			totalAnimTime = 1;
			
			% Time between frames
			dtFrame = totalAnimTime / max(nT - 1, 1);
	
			% find x limits and y limits for PDE plot
			xMin = min(self.domain.mesh.nodes(:,1));
			xMax = max(self.domain.mesh.nodes(:,1));
			yMin = min(self.domain.mesh.nodes(:,2));
			yMax = max(self.domain.mesh.nodes(:,2));
			heightRatio = (yMax - yMin) / (xMax - xMin);
			width = 500;
		
			% Create a new figure and fix its position
			%f = figure('Position',[100, 100, 2 * width, width * heightRatio]);

			% plot
			% NOTE: on 4/24/2025, I changed the order of the plots. Feel free to change back.
			% NOTE: on 4/24/2025, I changed this to have 3 subplots so I could manually add a frostbite plot for figures. remove this. 
			subplot(1,3,1);
			plot2(timestep);
			subplot(1,3,2);
			plot1(timestep);


			% adjust position
			f.Position = [100, 100, 2 * width, width * heightRatio];

			%{
			% Loop over timesteps
			for t = timesteps:timesteps

				% reset figure
				clf;             

				% plot current timestep
				subplot(1,2,1);
				plot1(t);
				subplot(1,2,2);
				plot2(t);

				% fix position
				f.Position = [100, 100, 2 * width, width * heightRatio];

				% delay
				pause(dtFrame);  
			end
			%}

		end


		% ANIMATION FUNCTIONS
		function animate(self,mode)
		% ANIMATE Sequentially plots the PDE + ODE solutions from the first
		% to the last timestep. The total animation duration is one second.
				
			if nargin == 1
				mode.plot1 = 'PDE';
				mode.plot2 = 'Temperatures';
			end

			% set first plot function
			if strcmp(mode.plot1,'PDE')
				plot1 = @(t) self.plotPDE(t);
			elseif strcmp(mode.plot1,'Frostbite')
				plot1 = @(t) self.plotFrostbite(t);
			elseif strcmp(mode.plot1,'BloodSupply')
				plot1 = @(t) self.plotBloodSupply(t);
			end

			% set second plot function
			if strcmp(mode.plot2,'ODE')
				plot2 = @(t) self.plotODE(t);
			elseif strcmp(mode.plot2,'Temperatures')
				plot2 = @(t) self.plotTemperatures(t);
			elseif strcmp(mode.plot2,'Average')
				plot2 = @(t) self.plotAverageSolutionValue(t);
			end

			% Number of timesteps
			nT = self.domain.time.N_t;
			
			% Total animation time (in seconds)
			totalAnimTime = 1;
			
			% Time between frames
			dtFrame = totalAnimTime / max(nT - 1, 1);
	
			% find x limits and y limits for PDE plot
			xMin = min(self.domain.mesh.nodes(:,1));
			xMax = max(self.domain.mesh.nodes(:,1));
			yMin = min(self.domain.mesh.nodes(:,2));
			yMax = max(self.domain.mesh.nodes(:,2));
			heightRatio = (yMax - yMin) / (xMax - xMin);
			width = 500;
		
			% Create a new figure and fix its position
			f = figure('Position',[100, 100, 2 * width, width * heightRatio]);

			% plot first time step
			subplot(1,2,1);
			plot1(1);
			subplot(1,2,2);
			plot2(1);

			% adjust position
			f.Position = [100, 100, 2 * width, width * heightRatio];
			pause();

			% Loop over timesteps
			for t = 1:nT

				% reset figure
				clf;             

				% plot current timestep
				subplot(1,2,1);
				plot1(t);
				subplot(1,2,2);
				plot2(t);

				% fix position
				f.Position = [100, 100, 2 * width, width * heightRatio];

				% delay
				pause(dtFrame);  
			end

		end
		
	end


	methods (Static)
		function c = redblue(m)
			%REDBLUE    Shades of red and blue color map
			%   REDBLUE(M), is an M-by-3 matrix that defines a colormap.
			%   The colors begin with bright blue, range through shades of
			%   blue to white, and then through shades of red to bright red.
			%   REDBLUE, by itself, is the same length as the current figure's
			%   colormap. If no figure exists, MATLAB creates one.
			%
			%   For example, to reset the colormap of the current figure:
			%
			%             colormap(redblue)
			%
			%   See also HSV, GRAY, HOT, BONE, COPPER, PINK, FLAG, 
			%   COLORMAP, RGBPLOT.
			%   Adam Auton, 9th October 2009
			if nargin < 1, m = size(get(gcf,'colormap'),1); end
			if (mod(m,2) == 0)
				% From [0 0 1] to [1 1 1], then [1 1 1] to [1 0 0];
				m1 = m*0.5;
				r = (0:m1-1)'/max(m1-1,1);
				g = r;
				r = [r; ones(m1,1)];
				g = [g; flipud(g)];
				b = flipud(r);
			else
				% From [0 0 1] to [1 1 1] to [1 0 0];
				m1 = floor(m*0.5);
				r = (0:m1-1)'/max(m1,1);
				g = r;
				r = [r; ones(m1+1,1)];
				g = [g; 1; flipud(g)];
				b = flipud(r);
			end
			c = [r g b]; 
		end
	end
end
