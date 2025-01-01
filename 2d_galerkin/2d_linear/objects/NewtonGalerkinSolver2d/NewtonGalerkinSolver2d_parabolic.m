classdef NewtonGalerkinSolver2d_parabolic < GalerkinSolver2d_parabolic & NewtonGalerkinSolver2d

    properties
        U
        iterHistory
    end

    methods
        function self = NewtonGalerkinSolver2d_parabolic(dom, auxfun)

            % Call superclass constructor
            self@GalerkinSolver2d_parabolic(dom, auxfun);
            self@NewtonGalerkinSolver2d();

        end

        function self = solveTimestep(self,S,b,FreeNodes)

			dirichlet = unique(self.domain.boundary.D_nodes);

            % define U_tilde (temporary)
            U_tilde = self.U;
            U_tilde(dirichlet) = 0;

            % Newton-Galerkin loop
            for iter = 1:100

                % Assembly
                self = self.assembleTensors;
                self = self.assembleVectors;
                self = self.assembleBCs;
                [DJ, J, S] = self.finalAssembly(self.U);

                %{
                % Solving one Newton step
                W = zeros(self.domain.mesh.nNodes,1);
                %W(unique(dirichlet)) = self.vectors.U_D(unique(dirichlet));
                W(FreeNodes) = DJ(FreeNodes,FreeNodes) \ J(FreeNodes);
                %W = DJ \ J;
                self.U(FreeNodes) = self.U(FreeNodes) - W(FreeNodes);
                %}

                % Solving one Newton step (temporary)
                W = zeros(self.domain.mesh.nNodes,1);
                W(FreeNodes) = DJ(FreeNodes,FreeNodes) \ J(FreeNodes);
                %W = DJ \ J;
                U_tilde = U_tilde - W;

                % check convegence
                if norm(W) < 10^(-10)
                    %fprintf(' %d iterations,',iter)
                    self.iterHistory(self.timestep) = iter;
                    break
                end
            end

            self.U = U_tilde + self.vectors.U_D;
            self.solution(:,self.timestep) = self.U;

            % store result
            %self.U = self.U + self.vectors.U_D;
            %self.solution(:,self.timestep) = self.U;
            %self.solution(:,self.timestep) = self.U + self.vectors.U_D;
        end

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

    end
end

        %{
        function self = solve(self)

            % Initialisation
            FreeNodes = self.domain.boundary.freeNodes;
            self = self.initializeProblem;

            %{
            % Update initial guess with Dirichlet BC values
            self = self.assembleBCs;
            self.U(unique(dirichlet)) = self.vectors.U_D(unique(dirichlet));
            %}

            %{
            % Time-stepping parameters
            t0 = 0; % Initial time
            tf = 1; % Final time
            dt = 0.01; % Time step size
            t = t0:dt:tf;
            %}

            % time-stepping loop
            for timestep = 1:self.domain.time.M_t

				self.timestep = timestep;
				self = self.initializeTimestep;
			    dirichlet = self.domain.boundary.D_nodes;
			    %self.U(unique(dirichlet)) = self.vectors.U_D(unique(dirichlet));

                % Newton-Galerkin loop
                for i = 1:100
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
                        fprintf(' %d iterations,',i)
                        break
                    end

                end

			    self = self.cleanup;

                % Update the solution for the next time step
                %self.U = self.U + dt * self.U;
            end
        end
        %}
