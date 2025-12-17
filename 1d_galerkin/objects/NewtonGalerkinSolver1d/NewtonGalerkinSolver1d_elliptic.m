classdef NewtonGalerkinSolver1d_elliptic < GalerkinSolver1d_elliptic & NewtonGalerkinSolver1d

	properties
	end

	methods
		function self = NewtonGalerkinSolver1d_elliptic(dom,auxfun)
			
			% call superclass constructor
			self@GalerkinSolver1d_elliptic(dom,auxfun);
			self@NewtonGalerkinSolver1d();

			self.U = ones(self.domain.mesh.nNodes,1);

		end

		function self = solve(self)

			%{
			% note: initial guess now set on instantiation
			%	This allows to change the initial guess before running solve.
			% Initial guess
			self.U = ones(self.domain.mesh.nNodes,1);
			%}

			% Initialization
			dirichlet = self.domain.boundary.D_nodes;
			FreeNodes = self.domain.boundary.freeNodes;
			self = self.assembleTensors;
			self = self.assembleVectors;

			% Update initial guess with Dirichlet BC values
			self = self.assembleBCs;
			self.U(unique(dirichlet)) = self.vectors.U_D(unique(dirichlet));

            % construct U_tilde
            U_tilde = self.U;
            U_tilde(dirichlet) = 0;

			% Newton-Galerkin iteration
			for i = 1:100
			
				% Assembly
				self = self.assembleTensors;
				self = self.assembleVectors;
				self = self.assembleBCs;
				[DJ,J] = self.finalAssembly(U_tilde);

				% Solving one Newton step
				W = zeros(self.domain.mesh.nNodes,1);
				W(FreeNodes) = DJ(FreeNodes,FreeNodes) \ J(FreeNodes);
                U_tilde = U_tilde - W;

                % add back dirichlet values to get current solution
				% note: I suspect this is necessary as the current solution will
				% be used to compute other tensors. Testing shows this gives the
				% best result in terms of convergence rate and accuracy.
				self.U = U_tilde + self.vectors.U_D;

				% check convegence
				if norm(W) < 10^(-10)
					%fprintf(' %d iterations,',i)
					break
				end

			end

			% store result
			self.solution = self.U;

		end


		function [DJ,J] = finalAssembly(self,U)

			% note: placeholder function. Actually handled by specific subclasses.
			% ...

		end

		function self = assembleTensors(self) 

			% call superclass method
			self = self.assembleTensors@GalerkinSolver1d_elliptic;

			% assemble additional tensors
			self.tensors.M_dr = self.assembleMassMatrix(self.coefficients.dr_du);
			self.tensors.M_dneu = self.computeNonlinearNeumannContribution;
			self.tensors.M_drob = self.computeNonlinearRobinContribution;

		end

		function self = assembleVectors(self) 

			% call superclass method
			self = self.assembleVectors@GalerkinSolver1d_elliptic;

			% assemble additional tensors
			self.vectors.M_r = self.computeVolumeForces(self.coefficients.r);

		end

		function self = assembleBCs(self)

			% note: as of Oct 11, 2024 periodic BCs don't work, so they are muted
			%self = self.computePeriodicBCs;
			self.vectors.U_D   = self.computeDirichletBCs;
			self.vectors.b_neu = self.computeNeumannBCs;
			%[self.tensors.M_rob,self.vectors.b_rob] = self.computeRobinBCs;
			%[self.tensors.M_rob,self.tensors.M_dyn,self.vectors.b_dyn] = self.computeDynamicBCs;

		end

	end
end


