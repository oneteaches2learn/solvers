classdef GalerkinSolver2d_heat < GalerkinSolver2d_parabolic

	methods
		function self = GalerkinSolver2d_heat(dom,cofs)
			
			% call superclass constructor
			self@GalerkinSolver2d_parabolic(dom,cofs);

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
