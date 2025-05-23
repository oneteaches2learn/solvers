classdef GalerkinSolver2d_heat < GalerkinSolver2d_parabolic

	methods
		function self = GalerkinSolver2d_heat(dom,auxfun)
			
			% call superclass constructor
			self@GalerkinSolver2d_parabolic(dom,auxfun);

			% solve
			self = self.solve();

		end

		function [S,b] = finalAssembly(self)

			% store variables
			tensors = self.tensors;
			vectors = self.vectors;

			% assemble LHS
			S = self.domain.time.dt * (tensors.A + tensors.M_rob + tensors.M_dyn_u) + ...
							tensors.M_p + tensors.M_dyn_du;

			% assemble RHS
			b = self.domain.time.dt * (vectors.b_vol - vectors.b_neu + vectors.b_rob + vectors.b_dyn) - ... 
					S * vectors.U_D + (tensors.M_p_prevTime + tensors.M_dyn_prevTime) * vectors.U_prevTime;

		end


	end
end
