classdef NewtonGalerkinSolver2d_heat < NewtonGalerkinSolver2d_parabolic

	properties
	end

	methods
		function self = NewtonGalerkinSolver2d_heat(dom,auxfun)
			
			% call superclass constructor
			self@NewtonGalerkinSolver2d_parabolic(dom,auxfun);

		end

		function [DJ,J] = finalAssembly(self,U)

            tensors = self.tensors;
            vectors = self.vectors;


			%{
			% copied from heat solver
			% assemble LHS
			S = self.domain.time.dt * (tensors.A + tensors.M_rob + tensors.M_dyn_u) + ...
							tensors.M_p + tensors.M_dyn_du;

			% assemble RHS
			b = self.domain.time.dt * (vectors.b_vol - vectors.b_neu + vectors.b_rob + vectors.b_dyn) - ... 
					S * vectors.U_D + (tensors.M_p_prevTime + tensors.M_dyn_prevTime) * vectors.U_prevTime;
			%}
			
			S = self.domain.time.dt * tensors.A + tensors.M_p;

			b = self.domain.time.dt * (vectors.b_vol) - ... 
					S * vectors.U_D + (tensors.M_p_prevTime) * vectors.U_prevTime;

			b = self.domain.time.dt * (vectors.b_vol) ... 
					+ (tensors.M_p_prevTime) * vectors.U_prevTime;

            J = (tensors.M_p + self.domain.time.dt * (tensors.A)) * self.U - b;
                
			DJ = tensors.M_p + self.domain.time.dt * tensors.A;
			



		end


	end
end
