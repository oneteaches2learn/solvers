classdef NewtonGalerkinSolver1d_poisson < NewtonGalerkinSolver1d_elliptic

	properties
	end

	methods
		function self = NewtonGalerkinSolver1d_poisson(dom,auxfun)
			
			% call superclass constructor
			self@NewtonGalerkinSolver1d_elliptic(dom,auxfun);

			% solve
			%self = self.solve;

		end

		% note: any mention of Robin BCs has been removed for now
		function [DJ,J] = finalAssembly(self,U_tilde)

			% store variables
            tensors = self.tensors;
            vectors = self.vectors;

			% assemble Linear Tensor
			%S = tensors.A + tensors.M_rob;
			S = tensors.A;

			% assemble nonlinear contributions to J
			b_nonlinear = self.vectors.M_r;

			% assemble Load Vector
			%b = (vectors.b_vol - vectors.b_neu + vectors.b_rob) - ... 
			%		S * vectors.U_D;
			b = (vectors.b_vol - vectors.b_neu) - ... 
					S * vectors.U_D;

			% assemble J
            J = S * U_tilde + b_nonlinear - b;
                
			% assemble DJ
			%DJ = tensors.A + tensors.M_rob + tensors.M_dr + tensors.M_dneu - tensors.M_drob;
			DJ = tensors.A + tensors.M_dr + tensors.M_dneu;


			%{
			% note: should self.tensors.M_rob be part of this? Investigate.
			% note: why should M_drob be subtracted? Investigate.
			DJ = self.tensors.A + self.tensors.M_dr + self.tensors.M_dneu + ...
						self.tensors.M_rob - self.tensors.M_drob;
			%J = (self.tensors.A + self.tensors.M_rob) * U + self.vectors.M_r - ...
			%				self.vectors.b_vol + self.vectors.b_neu - self.vectors.b_rob;

			J = (self.tensors.A + self.tensors.M_rob) * U + self.vectors.M_r - ...
							(self.vectors.b_vol + self.vectors.b_rob - self.vectors.b_neu);
			%}
		end

	end
end
