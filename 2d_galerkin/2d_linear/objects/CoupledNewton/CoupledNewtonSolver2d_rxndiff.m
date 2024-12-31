classdef CoupledNewtonSolver2d_rxndiff < NewtonGalerkinSolver2d_rxndiff

	properties
		ODE
		V_prev % value of V at previous time step
		V	   % value of V at current time and current iteration
	end

	methods
		function self = CoupledNewtonSolver2d_rxndiff(dom,auxfun,ODE)
			
			% call superclass constructor
			self@NewtonGalerkinSolver2d_rxndiff(dom,auxfun);

			% store additional data
			self.ODE = ODE;

			% solve
			self = self.solve;

		end

		function self = solve(self)

			% initialize problem
			FreeNodes = self.domain.boundary.freeNodes;
			self = self.initializeProblem;

			% initialize ODE (should be added to initializeTimestep eventually)
			self.V_prev = self.ODE.vInit;
			self.V = self.V_prev;
			self.ODE.solution(1) = self.V;

			for timestep = 1:self.domain.time.M_t

				% initialize timestep
				self.timestep = timestep;
				self = self.initializeTimestep;

				% initialize ODE (should be added to initializeTimestep eventually)
				self.V_prev = self.V;

				% assemble problem
				self  = self.assembleTensors;
				self  = self.assembleVectors;
				self  = self.assembleBCs;
				[S,b] = self.finalAssembly;

				% solve and store solution
				self = self.solveTimestep(S,b,FreeNodes);

				% break at equilibrium
				if self.equilibrium == 1, break; end

			end

			self = self.cleanup;

		end

        function self = solveTimestep(self,S,b,FreeNodes)

			dirichlet = self.domain.boundary.D_nodes;
   
            % Newton-Galerkin loop
            for iter = 1:100

                % Assembly
                self = self.assembleTensors;
                self = self.assembleVectors;
                self = self.assembleBCs;
                [DJ, J] = self.finalAssembly(self.U);

                % Dirichlet conditions
                W = zeros(self.domain.mesh.nNodes,1);
                W(unique(dirichlet)) = 0;
                
                % Solving one Newton step
                W(FreeNodes) = DJ(FreeNodes,FreeNodes) \ J(FreeNodes);
                self.U = self.U - W;

				% UPDATE ODE
				% store variables
				dt = self.ODE.dt;
				s  = self.ODE.s;
				R  = 1 / (1 + dt * s);
				g  = self.ODE.f;
				SU = self.domain.nodalQuadrature(self.U) / self.domain.domainArea;

				% solve current iteration
				self.V = R * (dt * g + dt * s * SU + self.V_prev); 

                % check convegence
                if norm(W) < 10^(-10)
                    %fprintf(' %d iterations,',iter)
                    self.iterHistory(self.timestep) = iter;
                    break
                end
            end

            % store result
            self.solution(:,self.timestep) = self.U + self.vectors.U_D;
			self.ODE.solution(self.timestep) = self.V;
        end

        function self = assembleVectors(self) 

            % call superclass method
            self = self.assembleVectors@NewtonGalerkinSolver2d_parabolic();

            % assemble additional tensors
            self.vectors.M_r = self.computeVolumeForces(self.coefficients.r);
			self.vectors.V   = 100 * self.V_prev;

        end

		%{
		function self = assembleTensors(self) 

			% call superclass method
			self = self.assembleTensors@NewtonGalerkinSolver2d_parabolic;

			% assemble additional tensors
			self.tensors.M_dr = self.assembleMassMatrix(self.coefficients.dr_du);

		end
		%}

		function [DJ,J] = finalAssembly(self,U)

			% store variables
            tensors = self.tensors;
            vectors = self.vectors;

			% assemble Linear Tensor
			S = self.domain.time.dt * (tensors.A + tensors.M_rob) + tensors.M_p;

			% assemble nonlinear contributions to J
			% RIGHT HERE IS THE CHANGE WHERE I AM INCORPORATING V
			b_nonlinear = self.domain.time.dt * (self.vectors.M_r - self.vectors.V);

			% assemble Load Vector
			b = self.domain.time.dt * (vectors.b_vol - vectors.b_neu + vectors.b_rob) - ... 
					S * vectors.U_D + (tensors.M_p_prevTime) * vectors.U_prevTime;

			% assemble J
            J = S * self.U + b_nonlinear - b;
                
			% assemble DJ
			DJ = tensors.M_p + self.domain.time.dt * (tensors.A + tensors.M_dr + tensors.M_dneu - tensors.M_drob);

		end
		%}

	end
end
