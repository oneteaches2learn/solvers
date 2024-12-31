classdef TimeStepping

	properties
		T
		M_t  % <~~ number of time steps
		N_t  % <~~ total number of time points
		dt
		equilibrium
	end

	properties (Dependent)
		tGrid
	end

	methods
		function self = TimeStepping(T,dt,eq)

			% set total time
			if nargin >= 1
				self.T = T;
			end

			% set dt and calculate M_t and N_t
			if nargin >= 2
				self.dt  = dt;
				%self.M_t = int32(T / dt);
				self.M_t = T / dt;
				self.N_t = self.M_t + 1;
			end

			% set equilibrium status
			if nargin == 3
				self.equilibrium = eq;
			else
				self.equilibrium = [];
			end

		end

		function self = setMesh(self,p,base)

			if nargin == 3
				% set dt from inputs
				self.dt = base^-p;
			end

			% set number of time steps
			self.M_t = self.T / self.dt;
			
			% check if M_t is an integer
			if mod(self.M_t,1) == 0
				...
			
			% else round N_t down to nearest integer
			else
				self.M_t = floor(self.N_t);
			end

			% set total number of time points
			self.N_t = self.M_t + 1;

		end


		% GETTERS
		function val = get.tGrid(self)

			val = linspace(0,self.T,self.N_t);
			
		end

	end
end
