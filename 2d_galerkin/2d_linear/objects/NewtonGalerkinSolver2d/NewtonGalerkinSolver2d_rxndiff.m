classdef NewtonGalerkinSolver2d_rxndiff < NewtonGalerkinSolver2d_parabolic

	properties
	end

	methods
		function self = NewtonGalerkinSolver2d_rxndiff(dom,auxfun)
			
			% call superclass constructor
			self@NewtonGalerkinSolver2d_parabolic(dom,auxfun);

			% solve
			% NOTE: Because the coupled solver is a subclass of this class,
			% calling solver here causes the coupled solver to attempt to solve
			% before the ODE is initialized. Eventually, you need to figure out
			% the class heirarchy structure that works for all solvers. Likely
			% this means that each class ending with a specific PDE name (e.g.
			% 'rxndiff') is a "leaf" and has no subclasses ever. An alternative
			% would be to take the "solve" method out entirely and call the
			% "solve" method from the outside. Or to implement an option to
			% solve or not solve in the constructor. I think the first option is
			% best. But I'm leaving this note for later consideration!  
			%self = self.solve;

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

		function [DJ,J] = finalAssembly(self,U_tilde)

			% store variables
            tensors = self.tensors;
            vectors = self.vectors;
			dt = self.domain.time.dt;

			% assemble Linear Tensor
			S = dt * (tensors.A + tensors.M_rob) + tensors.M_p;

			% assemble nonlinear contributions to J
			b_nonlinear = dt * self.vectors.M_r;

			% assemble Load Vector
			b = tensors.M_p_prevTime * vectors.U_prevTime + ... 
					dt * (vectors.b_vol - vectors.b_neu + vectors.b_rob) - ... 
					S * vectors.U_D;

			% assemble J
            J = S * U_tilde + b_nonlinear - b;

                
			% assemble DJ
			DJ = tensors.M_p + dt * (tensors.A + tensors.M_dr + tensors.M_dneu - tensors.M_drob);

		end

	end
end
