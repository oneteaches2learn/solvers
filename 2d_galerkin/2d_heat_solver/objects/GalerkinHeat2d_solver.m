classdef GalerkinHeat2d_solver < GalerkinParabolic2d_solver

	properties
	end

	methods
		function self = GalerkinHeat2d_solver(dom,time,cofs,uInit,f)
			
			% call superclass constructor
			self@GalerkinParabolic2d_solver(dom,time,cofs,uInit,f);

		end

		function tensors = initializeTensors(self)

			% create fields for tensor storage
			tensors.A   = [];
			tensors.M_p_prev = self.assembleMassMatrix(self.coefficients.p);
			tensors.M_p = [];
			tensors.E   = [];

			% check which tensors are time-varying
			tensors.timeVarying.A   = Coefficients.checkTimeVarying(self.coefficients.k);
			tensors.timeVarying.M_p = Coefficients.checkTimeVarying(self.coefficients.p);

		end

		function self = assembleTensors(self)

			% if first timestep, create tensors
			if self.isFirstTimeStep == 1
				self.tensors.A   = self.assembleStiffnessMatrix;
				self.tensors.M_p = self.assembleMassMatrix(self.coefficients.p);

			% else, update tensors as needed
			else

				% update A
				if self.tensors.timeVarying.A == 1
					self.tensors.A = self.assembleStiffnessMatrix;
				end

				% update M_p
				self.tensors.M_p_prev = self.tensors.M_p;
				if self.tensors.timeVarying.M_p == 1
					cof = self.coefficients.p;
					self.tensors.M_p = self.assembleMassMatrix(cof);
				end

			end

			% update Robin boundary tensor
			[self.tensors.E,temp] = self.computeRobinBCs;

		end

		function self = assembleVectors(self,t)

			% assemble vectors
			self.vectors.b_vol = self.computeVolumeForces(t);
			self.vectors.U_D   = self.computeDirichletBCs(t);
			self.vectors.b_neu = self.computeNeumannBCs(t);
			[temp,self.vectors.b_rob] = self.computeRobinBCs;

		end

		function [S,b] = finalAssembly(self,vectors,U_prev)

			% store variables
			tensors = self.tensors;
			vectors = self.vectors;

			% assemble LHS
			S = self.time.dt * (tensors.A + tensors.E) + tensors.M_p;

			% assemble RHS
			b = self.time.dt * (vectors.b_vol - vectors.b_neu + vectors.b_rob) - ... 
					S * vectors.U_D + tensors.M_p_prev * U_prev;

		end


	end
end
