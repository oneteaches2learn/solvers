classdef CoupledNewtonSolver2d_parabolic < NewtonGalerkinSolver2d_parabolic

    properties
        U
        iterHistory
    end

    methods
        function self = CoupledNewtonSolver2d_parabolic(dom, auxfun)

            % Call superclass constructor
            self@NewtonGalerkinSolver2d_parabolic(dom, auxfun);

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

                % check convegence
                if norm(W) < 10^(-10)
                    %fprintf(' %d iterations,',iter)
                    self.iterHistory(self.timestep) = iter;
                    break
                end
            end

            % store result
            self.solution(:,self.timestep) = self.U + self.vectors.U_D;
        end

        %{
        function self = initializeProblem(self)

            % call superclass method
            self = initializeProblem@GalerkinSolver2d_parabolic(self);

            % set initial guess
            self.U = self.solution(:,1);

        end

		function self = assembleTensors(self) 

			% call superclass method
			self = self.assembleTensors@GalerkinSolver2d_parabolic;

			% assemble additional tensors
			self.tensors.M_dneu = self.computeNonlinearNeumannContribution;
			self.tensors.M_drob = self.computeNonlinearRobinContribution;

		end
        %}

    end
end
