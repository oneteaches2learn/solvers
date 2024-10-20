classdef (Abstract) GalerkinSolver2d_parabolic < GalerkinSolver2d

	properties
		uInit
	end

	properties (Hidden)
		timestep
		t = 0
	end

	properties (Hidden,Dependent)
		isFirstTimeStep
	end

	methods
		function self = GalerkinSolver2d_parabolic(dom,auxfun)
			
			% call superclass constructor
			self@GalerkinSolver2d(dom,auxfun);

			if nargin == 2
				% store additional inputs
				self.uInit = auxfun.uInit;
				
				% calculate solution
				self = self.solve;
			end

		end

		function self = solve(self)

			% initialize problem
			FreeNodes = self.domain.boundary.freeNodes;
			self = self.initializeProblem;

			for timestep = 1:self.domain.time.M_t

				% initialize timestep
				self.timestep = timestep;
				self = self.initializeTimestep;

				% assemble problem
				self  = self.assembleTensors;
				self  = self.assembleVectors;
				self  = self.assembleBCs;
				[S,b] = self.finalAssembly;

				% solve and store solution
				v = sparse(self.domain.mesh.nNodes,1);
				v(FreeNodes) = S(FreeNodes,FreeNodes) \ b(FreeNodes);
				self.solution(:,self.timestep) = v + self.vectors.U_D;

				% break at equilibrium
				if self.equilibrium == 1, break; end

			end

			self = self.cleanup;

		end


		% INITIALIZATION FUNCTIONS
		function self = initializeProblem(self)

			% Record first timestep
			coords = self.domain.mesh.nodes;
			self.solution = zeros(self.domain.mesh.nNodes,self.domain.time.N_t);
			self.solution(:,1) = self.uInit(coords(:,1),coords(:,2));

			% Initialize tensors
			self = self.initializeTensors;
			self = self.initializeVectors;

		end

		function self = initializeTensors(self)

			% create fields for tensor storage
			tensors.A     = [];
			tensors.M_p   = [];
			tensors.M_rob = [];
			tensors.M_dyn = [];
			tensors.M_p_prevTime = self.assembleMassMatrix(self.coefficients.c);
			[temp,tensors.M_dyn_prevTime,temp] = self.computeDynamicBCs;

			% check which tensors are time-varying
			tensors.timeVarying.A   = Coefficients.isTimeVarying(self.coefficients.k);
			tensors.timeVarying.M_p = Coefficients.isTimeVarying(self.coefficients.c);

			% update tensors property
			self.tensors = tensors;

		end

		function self = initializeVectors(self)

			% create fields for vector storage
			self.vectors.b_vol = [];
			self.vectors.U_D   = [];
			self.vectors.b_neu = []; 
			self.vectors.b_rob = []; 
			self.vectors.b_dyn = []; 
			self.vectors.U_prevTime = []; 

			% check which vectors are time-varying
			self.vectors.timeVarying.b_vol = Coefficients.isTimeVarying(self.f);

		end

		function self = initializeTimestep(self)

			% store solution at previous timestep
			self.vectors.U_prevTime = self.solution(:,self.timestep);

			% update time stepping
			self.t = double(self.timestep) * double(self.domain.time.dt);
			self.timestep = self.timestep + 1;

		end


		% ASSEMBLY FUNCTIONS
		function self = assembleTensors(self)

			% if first timestep, create tensors
			if self.isFirstTimeStep == 1
				self.tensors.A   = self.assembleStiffnessMatrix;
				self.tensors.M_p = self.assembleMassMatrix(self.coefficients.c);

			% else, update tensors as needed
			else

				% update A
				if self.tensors.timeVarying.A == 1
					self.tensors.A = self.assembleStiffnessMatrix;
				end

				% update M_p
				if self.tensors.timeVarying.M_p == 1
					self.tensors.M_p_prevTime = self.tensors.M_p;
					cof = self.coefficients.c;
					self.tensors.M_p = self.assembleMassMatrix(cof);
				end

				% update M_dyn
				% note: have not yet implemented logic for time-varying dynamic BCs
				self.tensors.M_dyn_prevTime = self.tensors.M_dyn_du;

			end
		end

		function self = assembleVectors(self)

			% if first timestep, create vectors
			if self.isFirstTimeStep == 1
				self.vectors.b_vol = self.computeVolumeForces;

			% else, update vectors as needed
			else

				% update b_vol
				if self.vectors.timeVarying.b_vol == 1
					self.vectors.b_vol = self.computeVolumeForces;
				end
			end
		end

		function [S,b] = finalAssembly(self)

			% NOTE: placeholder function. Actually handled by specific subclasses.
			...

		end


		% UTILITY FUNCTIONS
		function result = equilibrium(self)

			% default value is zero
			result = 0;

			try % <~~~ if equilibrium is not set, then the code block below
				%	   produces an error. Using the "try" statement allows to
				%	   skip the following block of code in this case. 

				% check if solver should break
				if strcmp(self.domain.time.equilibrium.atEq,"break") == 1
					
					% get error between subsequent timesteps
					arg1 = self.solution(:,self.timestep);
					arg2 = self.solution(:,self.timestep-1);
					err = self.domain.L2err_threePointQuadrature_nodal(arg1,arg2);

					% if error < tolerance, break
					if err < self.domain.time.equilibrium.tolerance
						result = 1;
					end

					% issue `no equilibrium' warning
					if self.timestep == self.domain.time.N_t
						warn = " \n\tWARNING, Trial terminated without reaching equilibrium.\n ";
						fprintf(warn);
					end
					
				end
			end
		end

		function self = cleanup(self)

			% call superclass method
			self = self.cleanup@GalerkinSolver2d;

			% update stored parameters
			self.domain.time.N_t = self.timestep;
			self.domain.time.M_t = self.domain.time.N_t - 1;
			self.domain.time.T   = self.domain.time.M_t * self.domain.time.dt;
			self.solution = self.solution(:,1:self.domain.time.N_t);
			self.domain.time.equilibrium.t_eq = self.domain.time.M_t * self.domain.time.dt;

		end


		% GETTERS
		function val = get.isFirstTimeStep(self)

			val = (self.t == self.domain.time.dt);

		end


		% PLOTTING FUNCTIONS
		function p = plotMaxSolutionValue(self)

			% store variables
			maxVal = self.getMaxSolutionValue;
			timeGrid = self.domain.time.tGrid;

			% plot result
			p = plot(timeGrid,maxVal);

		end

		function p = plotMinSolutionValue(self)

			% store variables
			minVal = self.getMinSolutionValue;
			timeGrid = self.domain.time.tGrid;

			% plot result
			p = plot(timeGrid,minVal);

		end

		function p = plotAverageSolutionValue(self)

			% store variables
			avgVal = self.getAverageSolutionValue;
			timeGrid = self.domain.time.tGrid;

			% plot result
			p = plot(timeGrid,avgVal);

		end

		function animate(self)

			% store variables
			N_t = self.domain.time.N_t;

			% capture zLim bounds
			u_min = min(min(self.solution));
			u_max = max(max(self.solution));
			if u_max == u_min, u_max = u_min + 1; end

			% plot first time step
			self.plot(1);
			zlim([u_min u_max])
			pause();

			% animate remaining steps
			for i = 1:N_t
				self.plot(i);
				zlim([u_min u_max]);
				pause(1/N_t)
			end

		end

		function animatePatch(self)

			N_t = self.domain.time.N_t;

			% plot first time step
			self.plotPatch(1);
			pause();

			% animate remaining steps
			for i = 1:N_t
				self.plotPatch(i);
				pause(1/N_t)
			end

		end

		function animate_yline(self,NameValuePairs)

			arguments
				self
				NameValuePairs.plotInclusions = 'off'
			end

			N_t = self.domain.time.N_t;

			% plot first time step
			self.plot_yline(timestep=1,plotInclusions=NameValuePairs.plotInclusions);
			pause();

			% animate remaining steps
			for i = 1:N_t
				self.plot_yline(timestep=i,plotInclusions=NameValuePairs.plotInclusions);
				pause(1/N_t)
			end
		end



		% SOLUTION ANALYSIS
		function [t,timestep] = getConvergenceTime(self,tol)

			% default tolerance = 10^-6
			if nargin == 1, tol = 10^-6; end

			for n = 1:self.domain.time.M_t

				% store current and next timestep
				arg1 = self.solution(:,n);		
				arg2 = self.solution(:,n+1);		

				% compute error between subsequent timesteps
				err_n = arg1 - arg2;

				% compute L2 norm of error
				int = self.domain.L2norm_piecewiseLinear(err_n);

				% check against tolerance
				if int < tol
					timestep = n-1;
					t = self.domain.time.dt * timestep;
					return	
				end

			end

			% if loop completes, then the desired tolerance is never reached
			t = NaN;
			timestep = NaN;
			fprintf('WARNING: Tolerance not reached during time-series\n');
		end

	end
end
