classdef NewtonGalerkinSolver2d_poisson < NewtonGalerkinSolver2d_elliptic

	properties
	end

	methods
		function self = NewtonGalerkinSolver2d_poisson(dom,auxfun)
			
			% call superclass constructor
			self@NewtonGalerkinSolver2d_elliptic(dom,auxfun);

			% solve
			self = self.solve;

		end

		function [DJ,J] = finalAssembly(self,U)

			% note: should self.tensors.M_rob be part of this? Investigate.
			% note: why should M_drob be subtracted? Investigate.
			DJ = self.tensors.A + self.tensors.M_dr + self.tensors.M_dneu + ...
						self.tensors.M_rob - self.tensors.M_drob;
			%J = (self.tensors.A + self.tensors.M_rob) * U + self.vectors.M_r - ...
			%				self.vectors.b_vol + self.vectors.b_neu - self.vectors.b_rob;

			J = (self.tensors.A + self.tensors.M_rob) * U + self.vectors.M_r - ...
							(self.vectors.b_vol + self.vectors.b_rob - self.vectors.b_neu);
		end

	end
end
