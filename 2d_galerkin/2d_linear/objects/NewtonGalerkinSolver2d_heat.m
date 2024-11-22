classdef NewtonGalerkinSolver2d_heat < NewtonGalerkinSolver2d_parabolic

	properties
	end

	methods
		function self = NewtonGalerkinSolver2d_heat(dom,auxfun)
			
			% call superclass constructor
			self@NewtonGalerkinSolver2d_parabolic(dom,auxfun);

		end

		function [DJ,J] = finalAssembly(self,U)

			% store variables
            tensors = self.tensors;
            vectors = self.vectors;

			% assemble Tensor
			S = self.domain.time.dt * (tensors.A + tensors.M_rob) + tensors.M_p;

			% assemble Load Vector
			b = self.domain.time.dt * (vectors.b_vol - vectors.b_neu + vectors.b_rob) - ... 
					S * vectors.U_D + (tensors.M_p_prevTime) * vectors.U_prevTime;

			% assemble J
            J = S * self.U - b;
                
			% assemble DJ
			DJ = tensors.M_p + self.domain.time.dt * (tensors.A + tensors.M_dneu - tensors.M_drob);

		end

	end
end
