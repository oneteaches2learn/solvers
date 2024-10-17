classdef NewtonGalerkinSolver2d_poisson < NewtonGalerkinSolver2d_elliptic

	properties
	end

	methods
		function self = NewtonGalerkinSolver2d_poisson(dom,auxfun)
			
			% call superclass constructor
			self@NewtonGalerkinSolver2d_elliptic(dom,auxfun);

		end

		function [DJ,J] = finalAssembly(self,U)

			DJ = self.tensors.A + self.tensors.M_dr + self.tensors.M_dneu + ...
						self.tensors.M_rob - self.tensors.M_drob;
			J = (self.tensors.A + self.tensors.M_rob) * U + self.vectors.M_r - ...
							self.vectors.b_vol + self.vectors.b_neu - self.vectors.b_rob;

		end

	end
end
