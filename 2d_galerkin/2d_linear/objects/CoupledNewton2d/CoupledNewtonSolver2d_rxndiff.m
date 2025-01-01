classdef CoupledNewtonSolver2d_rxndiff < NewtonGalerkinSolver2d_rxndiff

	properties
		ODE
		V_prev % value of V at previous time step
		V	   % value of V at current time and current iteration
	end

	methods
		function self = CoupledNewtonSolver2d_rxndiff(dom,auxfun,ODE)
			
			% call superclass constructor
			self@NewtonGalerkinSolver2d_rxndiff(dom,auxfun);

			% store additional data
			self.ODE = ODE;
			self.ODE.tGrid = self.domain.time.tGrid;

			% solve
			self = self.solve;

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

			% temporary
			self.U = U_tilde;
            
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
				U_tilde(1)

				self.U = U_tilde;

                % add back dirichlet values to get current solution
				% note: I suspect this is necessary as the current solution will
				% be used to compute other tensors. But perhaps it is not
				% necessary, or even causes problems. Take a look.
				% self.U = U_tilde + self.vectors.U_D;

				%{
				% UPDATE ODE
				% note: If you were not timelagging the V value, you would update it here
				% store variables
				dt = self.ODE.dt;
				s  = self.ODE.s;
				R  = 1 / (1 + dt * s);
				g  = self.ODE.g;
				SU = self.domain.nodalQuadrature(self.U) / self.domain.domainArea;

				% solve current iteration
				self.V = R * (dt * g + dt * s * SU + self.V_prev); 
				%}

                % check convegence
                if norm(W) < 10^(-10)
                    fprintf(' %d iterations,\n',iter)
                    self.iterHistory(self.timestep) = iter;
                    break
                end
            end

			% (temporary) add back dirichlet values to get current solution
			self.U = U_tilde + self.vectors.U_D;

			% UPDATE ODE
			% Note: If you are time lagging the solutions, you would update the ODE here
			% store variables
			dt = self.ODE.dt;
			s  = self.ODE.s;
			R  = 1 / (1 + dt * s);
			g  = self.ODE.g;
			SU = self.domain.nodalQuadrature(self.U) / self.domain.domainArea;

			% solve current iteration
			self.V = R * (dt * g + dt * s * SU + self.V_prev); 

            % store resolved solution
            self.solution(:,self.timestep) = self.U;
			self.ODE.solution(self.timestep) = self.V;
        end

		function b = computeVolumeForces(self,f1,f2)

			% unpack variables
			nNodes    = self.domain.mesh.nNodes;
			elements3 = self.domain.mesh.elements(self.domain.mesh.effectiveElems, :);
			nElem3    = size(elements3,1);
			coords    = self.domain.mesh.nodes;

			% store function
			if nargin == 1,
				f1 = self.f;
				f2 = @(x1,x2,t)1; 
			elseif nargin == 2
				f2 = @(x1,x2,t)1;
			end

			% check coefficient variables
			[f1,t,U,V] = self.checkVariables(f1);
			f2 = self.checkVariables(f2);

			% initialize storage
			b = sparse(nNodes,1);

			% compute centroids of elements
			elementCoordsX = reshape(coords(elements3, 1), [nElem3, 3]);
			elementCoordsY = reshape(coords(elements3, 2), [nElem3, 3]);
			centroidsX = sum(elementCoordsX, 2) / 3;
			centroidsY = sum(elementCoordsY, 2) / 3;
			centroidsU = sum(U(elements3), 2) / 3;

			% compute volume forces
			areas = self.domain.mesh.areas(self.domain.mesh.effectiveElems);
			f1_vals = f1(centroidsX, centroidsY, t, centroidsU, V);
			f2_vals = f2(centroidsX, centroidsY, t, centroidsU, V);
			forces = f1_vals .* f2_vals / 3;
			volumeForces = areas .* forces;

			% accumulate forces into the global vector
			b = accumarray(elements3(:), repmat(volumeForces, 3, 1), [nNodes, 1], @sum, 0, true);

		end

		function [f,t,U,V] = checkVariables(self,f)

			% unpack variables
			nNodes    = self.domain.mesh.nNodes;

			% store function type data
			isCoupled     = num2str(Coefficients.isCoupled(f));
			isNonlinear   = num2str(Coefficients.isNonlinear(f));
			isTimeVarying = num2str(Coefficients.isTimeVarying(f));
			code = strcat(isCoupled,isNonlinear,isTimeVarying);
			code = str2num(code);

			% Update function signature and store variables and/or set dummy variables
			switch code
				% Spatially varying only
				case 000
					f = @(x1,x2,t,u,v)(f(x1,x2));
					t = 0;
					U = zeros(nNodes,1);
					V = 0;

				% Time varying
				case 001
					f = @(x1,x2,t,u,v)(f(x1,x2,t));
					t = self.t;
					U = zeros(nNodes,1);
					V = 0;
				
				% Nonlinear
				case 010
					f = @(x1,x2,t,u,v)(f(x1,x2,u));
					t = 0;
					U = self.U;
					V = 0;
				
				% Time varying and nonlinear
				case 011
					f = @(x1,x2,t,u,v)(f(x1,x2,t,u));
					t = self.t;
					U = self.U;
					V = 0;
				
				% Coupled
				case 100
					f = @(x1,x2,t,u,v)(f(x1,x2,v));
					t = 0;
					U = zeros(nNodes,1);
					%V = self.V;		% <~~~ use if V is not time-lagged
					V = self.V_prev;    % <~~~ use if V is time-lagged

				% Coupled and time varying
				case 101
					f = @(x1,x2,t,u,v)(f(x1,x2,t,v));
					t = self.t;
					U = zeros(nNodes,1);
					%V = self.V;		% <~~~ use if V is not time-lagged
					V = self.V_prev;    % <~~~ use if V is time-lagged

				% Coupled and nonlinear
				case 110
					f = @(x1,x2,t,u,v)(f(x1,x2,u,v));
					t = 0;
					U = self.U;
					%V = self.V;		% <~~~ use if V is not time-lagged
					V = self.V_prev;    % <~~~ use if V is time-lagged

				% Coupled, time varying, and nonlinear
				case 111
					f = @(x1,x2,t,u,v)(f(x1,x2,t,u,v));
					t = self.t;
					U = self.U;
					%V = self.V;		% <~~~ use if V is not time-lagged
					V = self.V_prev;    % <~~~ use if V is time-lagged
			end
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

			% Left panel: Use superclass's plot method to plot the PDE solution
			subplot(1, 2, 1);
			plot@GalerkinSolver2d(self,timestep);
			title(sprintf('u, t = %d',self.domain.time.tGrid(timestep)));

			% set zlim
			zMin = min(min(self.solution(:)),min(self.ODE.solution(:)));
			zMax = max(max(self.solution(:)),max(self.ODE.solution(:)));
			if zMin == zMax
				zMax = zMin + 1;
			end
			zlim([zMin, zMax]);

			% Right panel: Use the ODE's plot method to plot the ODE solution
			subplot(1, 2, 2);
			self.ODE.plot(timestep);
			title('v');

			% adjust figure position
			f = gcf;
			f.Position = [100, 100, 1000, 400];
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
			figure('Position',[100, 100, 1000, 400]);
			self.plot(1);
			pause();

			% Loop over timesteps
			for t = 1:nT
				clf;             % Clear the figure for redrawing
				self.plot(t);    % Call our custom plot method
				pause(dtFrame);  % Pause for the frame duration
			end
		end

	end
end
