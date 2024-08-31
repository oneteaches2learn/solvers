classdef GalerkinHeat2d_solver < GalerkinParabolic2d_solver

	methods
		function self = GalerkinHeat2d_solver(dom,cofs)
			
			% call superclass constructor
			self@GalerkinParabolic2d_solver(dom,cofs);

		end

		function [S,b] = finalAssembly(self)

			% store variables
			tensors = self.tensors;
			vectors = self.vectors;

			% assemble LHS
			S = self.domain.time.dt * (tensors.A + tensors.E) + tensors.M_p;

			% assemble RHS
			b = self.domain.time.dt * (vectors.b_vol - vectors.b_neu + vectors.b_rob) - ... 
					S * vectors.U_D + tensors.M_p_prevTime * vectors.U_prevTime;

		end


	end
end
