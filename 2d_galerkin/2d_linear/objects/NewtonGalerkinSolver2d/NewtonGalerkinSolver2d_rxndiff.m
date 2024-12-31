classdef NewtonGalerkinSolver2d_rxndiff < NewtonGalerkinSolver2d_parabolic

	properties
	end

	methods
		function self = NewtonGalerkinSolver2d_rxndiff(dom,auxfun)
			
			% call superclass constructor
			self@NewtonGalerkinSolver2d_parabolic(dom,auxfun);

		end

        function self = assembleVectors(self) 

            % call superclass method
            self = self.assembleVectors@NewtonGalerkinSolver2d_parabolic();

            % assemble additional tensors
            self.vectors.M_r = self.computeVolumeForces(self.coefficients.r);

        end

		function self = assembleTensors(self) 

			% call superclass method
			self = self.assembleTensors@NewtonGalerkinSolver2d_parabolic();

			% assemble additional tensors
			self.tensors.M_dr = self.assembleMassMatrix(self.coefficients.dr_du);

		end

		function [DJ,J] = finalAssembly(self,U)

			% store variables
            tensors = self.tensors;
            vectors = self.vectors;

			% assemble Linear Tensor
			S = self.domain.time.dt * (tensors.A + tensors.M_rob) + tensors.M_p;

			% assemble nonlinear contributions to J
			b_nonlinear = self.domain.time.dt * self.vectors.M_r;

			% assemble Load Vector
			b = self.domain.time.dt * (vectors.b_vol - vectors.b_neu + vectors.b_rob) - ... 
					S * vectors.U_D + (tensors.M_p_prevTime) * vectors.U_prevTime;

			% assemble J
            J = S * self.U + b_nonlinear - b;

                
			% assemble DJ
			DJ = tensors.M_p + self.domain.time.dt * (tensors.A + tensors.M_dr + tensors.M_dneu - tensors.M_drob);

		end

	end
end
