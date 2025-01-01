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

		function [DJ,J,S] = finalAssembly(self,U)

			% store variables
            tensors = self.tensors;
            vectors = self.vectors;
			dt = self.domain.time.dt;
			dirichlet = unique(self.domain.boundary.D_nodes);

			% temporary vectors
			U_tilde = self.U;
			U_tilde(dirichlet) = 0;
			U_tilde_prevTime = self.vectors.U_prevTime;
			U_tilde_prevTime(dirichlet) = 0;

			% assemble J (temporary)
			S = dt * tensors.A + tensors.M_p;
			b = tensors.M_p * U_tilde_prevTime + dt * vectors.b_vol - dt * S * vectors.U_D;
			J = S * U_tilde - b;

			% assemble DJ (temporary)
			DJ = tensors.M_p + self.domain.time.dt * (tensors.A + tensors.M_dneu - tensors.M_drob);

			%{
			% assemble Tensor
			S = dt * (tensors.A + tensors.M_rob) + tensors.M_p;

			% assemble Load Vector
			b = dt * (vectors.b_vol - vectors.b_neu + vectors.b_rob) - ... 
					dt * S * vectors.U_D + (tensors.M_p_prevTime) * vectors.U_prevTime;
			%b = self.domain.time.dt * (vectors.b_vol - vectors.b_neu + vectors.b_rob)  ... 
			%		+ (tensors.M_p_prevTime) * vectors.U_prevTime;

			% assemble J
            J = S * self.U - b;
                
			% assemble DJ
			DJ = tensors.M_p + self.domain.time.dt * (tensors.A + tensors.M_dneu - tensors.M_drob);
			%}

		end

	end
end
