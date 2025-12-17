classdef GalerkinSolver1d_poisson < GalerkinSolver1d_elliptic

    methods
        function self = GalerkinSolver1d_poisson(dom,auxfun)

            % call superclass constructor
            self@GalerkinSolver1d_elliptic(dom,auxfun);

        end
        
		function [S,b] = finalAssembly(self)

            % load variables
            tensors = self.tensors;
            vectors = self.vectors;

            % assemble LHS
			%S = tensors.A + tensors.M_r + tensors.M_rob;
			S = tensors.A + tensors.M_r;

            % assemble RHS
			%b = vectors.b_vol - vectors.b_neu + vectors.b_rob - S * vectors.U_D;
			b = vectors.b_vol - vectors.b_neu - S * vectors.U_D;

		end
    end

end