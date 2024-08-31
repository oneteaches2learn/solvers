classdef GalerkinRxndiff2d_solver < GalerkinParabolic2d_solver

	properties
	end

	methods
		function self = GalerkinRxndiff2d_solver(dom,cofs)
			
			% call superclass constructor
			self@GalerkinParabolic2d_solver(dom,cofs);

		end

		function self = initializeTensors(self,t)
			
			% call superclass method
			self = initializeTensors@GalerkinParabolic2d_solver(self);

			% create fields for tensor storage
			self.tensors.M_r = [];

			% check which tensors are time-varying
			self.tensors.timeVarying.M_r = Coefficients.isTimeVarying(self.coefficients.r);

		end

		function self = assembleTensors(self)

			% call superclass method
			self = assembleTensors@GalerkinParabolic2d_solver(self);

			% if first timestep, create tensors
			if self.isFirstTimeStep == 1
				self.tensors.M_r = self.assembleMassMatrix(self.coefficients.r);

			% else, update tensors as needed
			else
				% update M_r
				if self.tensors.timeVarying.M_r == 1
					cof = self.coefficients.r;
					self.tensors.M_r = self.assembleMassMatrix(cof);
				end
			end

		end

		function [S,b] = finalAssembly(self)

			% store variables
			tensors = self.tensors;
			vectors = self.vectors;

			% assemble LHS
			S = self.domain.time.dt * (tensors.A + tensors.M_r + tensors.E) + tensors.M_p;

			% assemble RHS
			b = self.domain.time.dt * (vectors.b_vol - vectors.b_neu + vectors.b_rob) - ...
				 				S * vectors.U_D + tensors.M_p_prevTime * vectors.U_prevTime;

		end

	end
end
