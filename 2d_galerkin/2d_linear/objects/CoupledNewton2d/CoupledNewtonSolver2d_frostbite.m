classdef CoupledNewtonSolver2d_frostbite < CoupledNewtonSolver2d_rxndiff
% This class is made to hold special data and functions specific to the
% frostbite model. The "solve" method is called by the superclass constructor.

	properties
	end

	methods
		function self = CoupledNewtonSolver2d_frostbite(dom,auxfun,ODE)

			% call superclass constructor
			self@CoupledNewtonSolver2d_rxndiff(dom,auxfun,ODE);

			% solve
			% ... called by superclass constructor

		end

		% PLOTTING FUNCTIONS
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
			t = self.domain.time.tGrid(timestep);
			Elements = self.domain.mesh.elements;
			Z = ones(size(X));

			% get solution values
			U = self.solution(:,timestep);
			V = self.ODE.solution(timestep);

			% get blood supply function
			r_activ = self.coefficients.r_activ;
			r_activ = self.checkVariables(r_activ);

			% compute U_data
			U_data = r_activ(X,Y,t,U,V) + zeros(size(X));
			U_data = U_data * nBuckets;
			U_data = round(U_data);

			% plot
			h = trisurf(Elements,X,Y,Z);
			view(2)
			str = sprintf('%% Blood Supply, t = %.2f',self.domain.time.tGrid(timestep));
			title(str);

			% convert U to length(U) x 3 matrix
			%U_data = reshape(U_normalized, [length(U_normalized), 1, 1]);
			C_data = self.redblue(nBuckets+1);

			% replace U_normalized entries with corresponding RGB values
			C = C_data(U_data+1,:);
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

		function animateBloodSupply(self)

			mode = struct('plot1','BloodSupply','plot2','Temperatures');
			self.animate(mode);
		
		end

		function plotFrostbite(self,timestep)
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

		function animateFrostbite(self)

			mode = struct('plot1','Frostbite','plot2','Temperatures');
			self.animate(mode);

		end

		function plotFrostbitePercent(self)
		% PLOTFROSTBITE Plots the mesh cells for the PDE solution at a given
		% timestep. Cells where U_centroids <= 0 are colored black; all others white.

			% get mesh data
			X = self.domain.mesh.nodes(:,1);
			Y = self.domain.mesh.nodes(:,2);
			Elements = self.domain.mesh.elements;
			Z = ones(size(X));

			for timestep = 1:self.domain.time.N_t
				% get U values at the centroids
				U = self.solution(:,timestep);
				U_centroids = self.domain.nodes2centroids(U);
				U_frostbite = find(U_centroids <= 0);

				% get area of frostbitten cells
				area_frostbite = sum(self.domain.mesh.areas(U_frostbite));
				frostbite_percent(timestep) = area_frostbite / sum(self.domain.mesh.areas);
			end

			t = self.domain.time.tGrid;

			% plot
			plot(t,frostbite_percent,'LineWidth',4);
			grid on
			title('Frostbite Area Percent','Interpreter','latex');
			xlabel('$t$','Interpreter','latex');
			ylabel('\%','Interpreter','latex');

		end

		function plot_sCoefficient(self)

			% store data
			tGrid = self.domain.time.tGrid;
			u = self.getAverageSolutionValue';
			v = self.ODE.solution;

			% compute and plot s values
			s = self.ODE.s(u,v) + zeros(size(tGrid));
			plot(tGrid,s,'LineWidth',4);
			grid on;

			% adjust ylim
			yMin = min(min(s));
			yMax = max(max(s));
			ylim([yMin - 0.1 * (yMax - yMin), yMax + 0.1 * (yMax - yMin)]);

			% add titles and labels
			xlabel('$t$','Interpreter','latex');
			ylabel('$s$','Interpreter','latex');
			title('$s(t)$','Interpreter','latex');

			% adjust figure position
			f = gcf;
			f.Position = [100, 100, 500, 600];

		end

	end
end
