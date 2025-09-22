classdef GalerkinSolver2d_poisson < GalerkinSolver2d_elliptic

    methods
        function self = GalerkinSolver2d_poisson(dom,auxfun)

            % call superclass constructor
            self@GalerkinSolver2d_elliptic(dom,auxfun);

        end
        
		function [S,b] = finalAssembly(self)

            % load variables
            tensors = self.tensors;
            vectors = self.vectors;

            % assemble LHS
			S = tensors.A + tensors.M_r + tensors.M_rob;

            % assemble RHS
			b = vectors.b_vol - vectors.b_neu + vectors.b_rob - S * vectors.U_D;

		end
    end

end