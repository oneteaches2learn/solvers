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

		function self = solve(self)

			% initialize problem
			self = self.initializeProblem;

			% initialize ODE (should be added to initializeTimestep eventually)
			self.V_prev = self.ODE.vInit;
			self.V = self.V_prev;
			self.ODE.solution = zeros(1,self.domain.time.N_t);
			self.ODE.solution(1) = self.V;

			for timestep = 1:self.domain.time.M_t

				% initialize timestep
				self.timestep = timestep;
				self = self.initializeTimestep;

				% initialize ODE (should be added to initializeTimestep eventually)
				self.V_prev = self.V;

				% solve and store solution
				self = self.solveTimestep();

				% break at equilibrium
				if self.equilibrium == 1, break; end

			end

			self = self.cleanup;

		end


        function self = solveTimestep(self,S,b,FreeNodes)

            % store dirichlet nodes and free nodes
			FreeNodes = self.domain.boundary.freeNodes;
			dirichlet = unique(self.domain.boundary.D_nodes);

            % construct U_tilde
            U_tilde = self.U;
            U_tilde(dirichlet) = 0;

            % Newton-Galerkin loop
            for iter = 1:10

                % Assembly
                self = self.assembleTensors;
                self = self.assembleVectors;
                self = self.assembleBCs;
                [DJ, J] = self.finalAssembly(U_tilde);

                % Solving one Newton step
                W = zeros(self.domain.mesh.nNodes,1);
                W(FreeNodes) = DJ(FreeNodes,FreeNodes) \ J(FreeNodes);
                U_tilde = U_tilde - W;

                % add back dirichlet values to get current solution
				% note: I suspect this is necessary as the current solution will
				% be used to compute other tensors. But perhaps it is not
				% necessary, or even causes problems. Take a look.
				self.U = U_tilde + self.vectors.U_D;

				% UPDATE ODE
				% note: If you were not timelagging the V value, you would update it here
				% store variables
				SU = self.domain.nodalQuadrature(self.U) / self.domain.domainArea;

				% remaining coefficients
				dt = self.ODE.dt;
				s  = self.ODE.s(SU,self.V);
				R  = 1 / (1 + dt * s);
				g  = self.ODE.g;

				% solve current iteration
				self.V = R * (dt * g + dt * s * SU + self.V_prev); 

				% apply resolvent
				vUpper = 37;
				self.V = self.ODE.resolvent(self.V, dt, s, NaN, vUpper);

                % check convegence
                if norm(W) < 10^(-10)
                    fprintf(' timestep: %d, converged in %d iterations\n',self.timestep,iter)
                    self.iterHistory(self.timestep) = iter;
                    break
                end
            end

			% (temporary) trying to figure out when v starts going back up?
			self.s_history(self.timestep) = s * (self.V - SU) / (1 + dt);

			%{
			% UPDATE ODE
			% Note: If you are time lagging the solutions, you would update the ODE here
			% store variables
			% compute S(u)
			SU = self.domain.nodalQuadrature(self.U) / self.domain.domainArea;

			%{
			% compute weighted S(u)
			U = self.U;
			U_centroids = self.domain.nodes2centroids(U);
			U_marked = find(U_centroids < 10);
			U_nan = U_centroids;
			U_nan(U_marked) = NaN;
			areas = self.domain.mesh.areas;
			areas_nan = areas;
			areas_nan(U_marked) = NaN;
			area_unmarked = sum(areas_nan);
			SU = self.domain.nodalQuadrature(U_nan) / area_unmarked;
			%}

			% remaining coefficients
			dt = self.ODE.dt;
			s  = self.ODE.s(SU,self.V);
			R  = 1 / (1 + dt * s);
			g  = self.ODE.g;

			% solve current iteration
			self.V = R * (dt * g + dt * s * SU + self.V_prev); 

			% apply resolvent
			vUpper = 37;
			self.V = self.ODE.resolvent(self.V, dt, s, NaN, vUpper);
			%}

            % store resolved solution
            self.solution(:,self.timestep) = self.U;
			self.ODE.solution(self.timestep) = self.V;
        end


		% PLOTTING FUNCTIONS
		function plot(self,timestep)
			% PLOT Plots the PDE solution on the left and the ODE solution on the right.
			% If an argument is passed, use that as the timestep; 
			% otherwise, use the last timestep from self.domain.time.tGrid.

			% default to last timestep
			if nargin == 1
				timestep = self.domain.time.N_t;
			end

			% calculate average hand temp
			SU = self.domain.nodalQuadrature(self.solution()) / self.domain.domainArea;

			% LEFT PANEL: Use superclass's plot method to plot the PDE solution
			%subplot(1, 2, 1);
			plot@GalerkinSolver2d(self,timestep);
			title(sprintf('u, t = %.2f',self.domain.time.tGrid(timestep)));
			view(2);

			% set zlim
			zMin = min([min(self.solution(:)),min(self.ODE.solution(:)),min(SU)]);
			zMax = max([max(self.solution(:)),max(self.ODE.solution(:)),max(SU)]);
			if zMin == zMax
				zMax = zMin + 1;
			end
			zlim([zMin, zMax]);

			% set xlim
			xMin = min(self.domain.mesh.nodes(:,1));
			xMax = max(self.domain.mesh.nodes(:,1));
			xlim([xMin, xMax]);

			% set ylim
			yMin = min(self.domain.mesh.nodes(:,2));
			yMax = max(self.domain.mesh.nodes(:,2));
			ylim([yMin, yMax]);


			% set colorbar
			colorbar;
			clim([zMin, zMax]);

			%{
			% RIGHT PANEL: Use the ODE's plot method to plot the ODE solution
			subplot(1, 2, 2);
			hold on
			h = self.ODE.plot(timestep);
			plot(self.ODE.tGrid(1:timestep), SU(1:timestep), 'LineWidth', 2);
			hold off

			% set zlim
			yMin = min([min(self.ODE.solution(:)),min(SU)]);
			yMax = max([max(self.ODE.solution(:)),max(SU)]);
			if yMin == yMax
				yMax = yMin + 1;
			end
			ylim([yMin - 0.1*yMin, yMax + 0.1*yMax]);

			title('v');

			% adjust figure position
			f = gcf;
			f.Position = [100, 100, 1000, 400];
			%}
			f = gcf;
			f.Position = [100, 100, 500, 600];
		end

		function h = plotTemperatures(self)

			% calculate average hand temp
			SU = self.domain.nodalQuadrature(self.solution()) / self.domain.domainArea;

			% plot
			hold on
			h = self.ODE.plot();
			t = self.ODE.tGrid();
            indSkip = ceil(length(t) / 10);
			plot(t, SU, '-s', 'LineWidth', 4, 'MarkerSize', 10, 'MarkerIndices', [1:indSkip:length(t)]);
			hold off

			% set title, legend, and labels
			title('Temperature vs. Time');
			xlabel('Time');
			ylabel('Temperature');
			legend('$v$','$\langle u \rangle_{\Omega}$','interpreter','latex','Location','best');

			% set zlim
			yMin = min([min(self.ODE.solution(:)),min(SU)]);
			yMax = max([max(self.ODE.solution(:)),max(SU)]);
			if yMin == yMax
				yMax = yMin + 1;
			end
			ylim([yMin - 0.1*yMin, yMax + 0.1*yMax]);


			% adjust figure position
			f = gcf;
			f.Position = [100, 100, 500, 500];

		end


		function h = plotFrostbite(self,timestep)
		% PLOTFROSTBITE Plots the mesh cells for the PDE solution at a given
		% timestep. Cells where U_centroids <= 0 are colored black; all others white.

			% Default to the last timestep
			if nargin == 1
				timestep = self.domain.time.N_t;
			end

			% get mesh data
			X = self.domain.mesh.nodes(:,1);
			Y = self.domain.mesh.nodes(:,2);
			Elements = self.domain.mesh.elements;
			Z = ones(size(X));

			% get U values at the centroids
			U = self.solution(:,timestep);
			U_centroids = self.domain.nodes2centroids(U);
			U_frostbite = find(U_centroids <= 0);

			% get nodes of frostbite cells
			U_nodes = Elements(U_frostbite,:);
			U_nodes = unique(U_nodes(:));
			
			% plot
			hold on 
			h1 = trisurf(Elements,X,Y,Z);
			h2 = trisurf(Elements,X,Y,Z);
			title(sprintf('Frostbitten Cells, t = %.2f',self.domain.time.tGrid(timestep)));
			view(2)
			hold off

			% create legend that indicates that frostbitten cells are black
			C1 = zeros([size(Z) 3]);	
			h1.CData = C1;
			legend(h1,'Frostbitten Cells','Location','northeast');
			k = gcf().CurrentAxes.Legend;
			k.AutoUpdate = 'off';

			% color the frostbite cells black
			C2 = ones([size(Z) 3]);	
			C2(U_nodes,:,:) = zeros([size(U_nodes) 3]);
			h2.CData = C2;

			% set xlim
			xMin = min(self.domain.mesh.nodes(:,1));
			xMax = max(self.domain.mesh.nodes(:,1));
			xlim([xMin, xMax]);

			% set ylim
			yMin = min(self.domain.mesh.nodes(:,2));
			yMax = max(self.domain.mesh.nodes(:,2));
			ylim([yMin, yMax]);

			% adjust figure position
			f = gcf;
			f.Position = [100, 100, 500, 600];

		end

		function plotBloodSupply(self,timestep)

			% (temporary) set number of buckets
			nBuckets = 10;

			% Default to the last timestep
			if nargin == 1
				timestep = self.domain.time.N_t;
			end

			% get mesh data
			X = self.domain.mesh.nodes(:,1);
			Y = self.domain.mesh.nodes(:,2);
			Elements = self.domain.mesh.elements;
			Z = ones(size(X));

			% get U values at the centroids
			U = self.solution(:,timestep);

			% clip U values above 32 and below 10
			U = max(U,10);
			U = min(U,32);

			% normalize U values
			U_normalized = (U - 10) / (32 - 10);
			U_normalized = round(U_normalized * nBuckets) / nBuckets;
			U_normalized = U_normalized * nBuckets;

			% plot
			h = trisurf(Elements,X,Y,Z);
			view(2)
			str = sprintf('%% Blood Supply, t = %.2f',self.domain.time.tGrid(timestep));
			title(str);

			% convert U to length(U) x 3 matrix
			U_data = reshape(U_normalized, [length(U_normalized), 1, 1]);
			C_data = self.redblue(nBuckets+1);

			% replace U_normalized entries with corresponding RGB values
			C = C_data(U_normalized+1,:);
			C = reshape(C, [length(C), 1, 3]);
			h.CData = C;

			% create colorbar
			colormap(C_data)
			colorbar;

			% set xlim
			xMin = min(self.domain.mesh.nodes(:,1));
			xMax = max(self.domain.mesh.nodes(:,1));
			xlim([xMin, xMax]);

			% set ylim
			yMin = min(self.domain.mesh.nodes(:,2));
			yMax = max(self.domain.mesh.nodes(:,2));
			ylim([yMin, yMax]);

			% adjust figure position
			f = gcf;
			f.Position = [100, 100, 500, 600];
			
		end

		function animate(self)
		% ANIMATE Sequentially plots the PDE + ODE solutions from the first
		% to the last timestep. The total animation duration is one second.
			
			% Number of timesteps
			nT = self.domain.time.N_t;
			
			% Total animation time (in seconds)
			totalAnimTime = 1;
			
			% Time between frames
			dtFrame = totalAnimTime / max(nT - 1, 1);
			
			% Create a new figure and fix its position
			figure('Position',[100, 100, 1000, 600]);
			self.plot(1);
			pause();

			% Loop over timesteps
			for t = 1:nT
				clf;             % Clear the figure for redrawing
				self.plot(t);    % Call our custom plot method
				pause(dtFrame);  % Pause for the frame duration
			end
		end

		function animateFrostbite(self)
			
			% Number of timesteps
			nT = self.domain.time.N_t;
			
			% Total animation time (in seconds)
			totalAnimTime = 1;
			
			% Time between frames
			dtFrame = totalAnimTime / max(nT - 1, 1);
			
			% Create a new figure and fix its position
			figure('Position',[100, 100, 500, 600]);
			self.plotFrostbite(1);
			pause();

			% Loop over timesteps
			for t = 1:nT
				clf;             
				self.plotFrostbite(t);
				pause(dtFrame);  
			end
		end

		function animateBloodSupply(self)
			
			% Number of timesteps
			nT = self.domain.time.N_t;
			
			% Total animation time (in seconds)
			totalAnimTime = 1;
			
			% Time between frames
			dtFrame = totalAnimTime / max(nT - 1, 1);
			
			% Create a new figure and fix its position
			figure('Position',[100, 100, 500, 600]);
			self.plotBloodSupply(1);
			pause();

			% Loop over timesteps
			for t = 1:nT
				clf;             
				self.plotBloodSupply(t);    
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
