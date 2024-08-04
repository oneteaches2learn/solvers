classdef GalerkinHeat2d_solver < GalerkinParabolic2d_solver

	properties
	end

	methods
		function self = GalerkinHeat2d_solver(dom,time,cofs,uInit,f)
			
			% call superclass constructor
			self@GalerkinParabolic2d_solver(dom,time,cofs,uInit,f);

		end

		function [S,b] = finalAssembly(self)

			% store variables
			tensors = self.tensors;
			vectors = self.vectors;

			% assemble LHS
			S = self.time.dt * (tensors.A + tensors.E) + tensors.M_p;

			% assemble RHS
			b = self.time.dt * (vectors.b_vol - vectors.b_neu + vectors.b_rob) - ... 
					S * vectors.U_D + tensors.M_p_prevTime * vectors.U_prevTime;

		end


	end
end
