classdef NewtonGalerkinSolver2d_poisson < NewtonGalerkinSolver2d_elliptic

	properties
	end

	methods
		function self = NewtonGalerkinSolver2d_poisson(dom,auxfun)
			
			% call superclass constructor
			self@NewtonGalerkinSolver2d_elliptic(dom,auxfun);

		end

		function [DJ,J] = finalAssembly(self,U)

			% Old way
			DJ = self.tensors.A + self.tensors.M_dr + self.tensors.M_dneu + ...
						self.tensors.M_rob - self.tensors.M_drob;
			J = (self.tensors.A + self.tensors.M_rob) * U + self.vectors.M_r - ...
							self.vectors.b_vol + self.vectors.b_neu - self.vectors.b_rob;

			%{
			DJ = self.finalStiffnessMatrixDerivative;
			J =  self.finalStiffnessMatrix * U - self.finalLoadVector;
			%}

		end

		function S = finalStiffnessMatrix(self)

			% load variables
			tensors = self.tensors;

			% assemble LHS
			S = tensors.M_p + tensors.A + tensors.M_rob;

		end

		function dS = finalStiffnessMatrixDerivative(self)

			% load variables
			tensors = self.tensors;

			% assemble LHS derivative
			dS = tensors.A + tensors.M_dr + tensors.M_dneu + tensors.M_rob - tensors.M_drob;

		end

		function b = finalLoadVector(self)

			% load variables
			vectors = self.vectors;

			% assemble RHS
			b = vectors.b_vol - vectors.b_neu + vectors.b_rob - vectors.M_r;
			% note: I don't know what M_r is or where it comes from. Investigate.

		end

	end
end
