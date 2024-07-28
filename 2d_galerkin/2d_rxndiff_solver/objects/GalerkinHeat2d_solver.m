classdef GalerkinHeat2d_solver < GalerkinParabolic2d_solver

	properties
	end

	methods
		function self = GalerkinHeat2d_solver(dom,time,cofs,uInit,f)
			
			% call superclass constructor
			self@GalerkinParabolic2d_solver(dom,time,cofs,uInit,f);

		end

		function tensors = assembleTensors(self,t)

			% assemble tensors
			tensors.A   = self.assembleStiffnessMatrix(t);
			tensors.M_r = self.assembleMassMatrix(self.coefficients.r,t);
			tensors.M_p = self.assembleMassMatrix(self.coefficients.p,t);
			[tensors.E,temp] = self.computeRobinBCs(t);

		end

		function vectors = assembleVectors(self,t)

			% assemble vectors
			vectors.b_vol = self.computeVolumeForces(t);
			vectors.U_D   = self.computeDirichletBCs(t);
			vectors.b_neu = self.computeNeumannBCs(t);
			[temp,vectors.b_rob] = self.computeRobinBCs(t);

		end

		function [S,b] = finalAssembly(self,tensors,vectors,U_prev)

			% assemble LHS
			S = self.time.dt * (tensors.A + tensors.M_r + tensors.E) + tensors.M_p;

			% assemble RHS
			b = self.time.dt * (vectors.b_vol - vectors.b_neu + vectors.b_rob) - ... 
					S * vectors.U_D + tensors.M_p * U_prev;

		end


	end
end
