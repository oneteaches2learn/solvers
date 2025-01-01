classdef NewtonGalerkinSolver2d_heat < NewtonGalerkinSolver2d_parabolic

	properties
	end

	methods
		function self = NewtonGalerkinSolver2d_heat(dom,auxfun)
			
			% call superclass constructor
			self@NewtonGalerkinSolver2d_parabolic(dom,auxfun);

			% solve
			self = self.solve;

		end

		function [DJ,J,S] = finalAssembly(self,U_tilde)

			% store variables
            tensors = self.tensors;
            vectors = self.vectors;
			dt = self.domain.time.dt;

			% assemble linear tensor
			S = dt * (tensors.A + tensors.M_rob) + tensors.M_p;

			% assemble Load Vector
			b = tensors.M_p_prevTime * vectors.U_prevTime + ... 
					dt * (vectors.b_vol - vectors.b_neu + vectors.b_rob) - ... 
					S * vectors.U_D;

			% assemble J
			J = S * U_tilde - b;

			% assemble DJ
			DJ = tensors.M_p + dt * (tensors.A + tensors.M_dneu - tensors.M_drob);

		end

	end
end
