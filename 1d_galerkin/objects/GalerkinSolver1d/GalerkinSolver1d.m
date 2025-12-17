classdef GalerkinSolver1d

    properties
        domain
        coefficients
        f
        solution
    end

    properties (Hidden)
        tensors
        vectors
        t = 0;    % current time
        U   % current solution
        V   % current coupled solution
        U_D % Dirichlet BC values
        U_N % Neumann BC values
        U_R % Robin BC values
    end

    methods 
        % CONSTRUCTOR
        function self = GalerkinSolver1d(dom,auxfun)

            if nargin == 2

                % store inputs
                self.domain = dom;
                self.coefficients = auxfun.cofs;
                self.f = auxfun.f;

            end

            % calculate solution
            % ...

        end

        % ASSEMBLY FUNCTIONS
        function A = assembleStiffnessMatrix(self)

            A_vec = self.localStiffnessMatrix_vec;
            A = self.sumLocalMatrices(A_vec);

        end

        function M = assembleMassMatrix(self,c)
            
            M_vec = self.localMassMatrix_vec(c);
            M = self.sumLocalMatrices(M_vec);

        end

        function A_vec = localStiffnessMatrix_vec(self,k)
        % LOCALSTIFFNESSMATRIX_VEC performs vectorized computation of 1D P1
        %   local stiffness matrices.  
        %
        % USAGE:
        %   A_vec = localStiffnessMatrix_vec(self,k)
        % INPUTS:
        %   k : (optional) function handle, diffusion coefficient
        % OUTPUTS:
        %   A_vec : 2 x 2 x nElems array
        %       A_vec(:,:,e) is the 2x2 local stiffness matrix for element e.

            % input handling
            if nargin < 2, k = self.coefficients.k; end

            % unpack variables
            elems  = self.domain.mesh.elements;  
            mids   = self.domain.mesh.midpoints;
            nElems = size(elems, 1);
            areas  = self.domain.mesh.areas;

            % check coefficient variables
            [k, t, U, V] = self.checkVariables(k);

            % compute k on element midpoints
            K = k(mids, t, U, V);     
            if length(K) == 1
                K = K * ones(nElems, 1); 
            end

            % get nElems x 1 vector of scale factors
            scale = K ./ areas;

            % reshape to 1 x 1 x nElems stack
            scale = reshape(scale, 1, 1, []);   % 1x1x nElems

            % store base stiffness pattern for 1D P1 method
            base = [ 1, -1; -1,  1 ]; 

            % Final result: A_vec(i,j,e) = scale(e) * base(i,j)
            A_vec = base .* scale;      

        end

        function M_vec = localMassMatrix_vec(self, c)
        % LOCALMASSMATRIX_VEC implements vectorized computation of 1D P1 local
        %   mass matrices.  
        % USAGE:
        %   M_vec = localMassMatrix_vec(self,c)
        % INPUTS:
        %   c : mass coefficient (function_handle)
        % OUTPUTS:
        %   M_vec : 2 x 2 x nElems array
        %       M_vec(:,:,e) is the 2x2 local mass matrix for element

            % unpack variables
            coords = self.domain.mesh.nodes(:);   
            mids   = self.domain.mesh.midpoints;
            elems  = self.domain.mesh.elements;
            nElems = size(elems, 1);
            areas  = self.domain.mesh.areas;

            % check coefficient variables
            [c, t, U, V] = self.checkVariables(c);
            
            % interpolate U on element midpoints
            U_centroids = mean(U(elems), 2);   % nElems x 1

            % evaluate c at midpoints
            C = c(mids, t, U_centroids, V);
            if numel(C) == 1
                C = C * ones(nElems, 1); 
            end

            % combine coefficient with element lengths (i.e. areas)
            Ch = C .* areas;

            % compute local mass matrices
            M_loc = self.localMassMatrix();    % 2 x 2
            Ch_vec = permute(Ch, [3, 2, 1]);   % 1 x 1 x nElems
            M_vec = pagemtimes(M_loc, Ch_vec); % 2 x 2 x nElems

        end

        function M = localMassMatrix(self)
            % localMassMatrix (1D reference)
            %
            % Consistent mass matrix for P1 on reference interval [0,1]:
            % ∫_0^1 [φ1 φ2]^T [φ1 φ2] dx = (1/6)[2 1; 1 2]

            M = [2, 1;
                1, 2] / 6;
        end

        function A = sumLocalMatrices(self,A_vec)
        % sumLocalMatrices Assemble global stiffness/mass matrix from
        %   vectorized local element matrices.
        % USAGE:
        %   A = sumLocalMatrices(self,A_vec)
        % INPUTS:
        %   A_vec : 2 x 2 x nElems array of local element matrices
        % OUTPUTS:
        %   A : sparse global matrix

            % unpack mesh connectivity variables
            coords = self.domain.mesh.nodes(:); 
            nNodes = numel(coords);
            elems  = self.domain.mesh.elements;

            % extract local entries as column vectors, 
            %   i.e. do A_vec(i,j,e) -> Aij(e)
            A11 = reshape(A_vec(1,1,:), [], 1);   % nElems x 1
            A12 = reshape(A_vec(1,2,:), [], 1);
            A21 = reshape(A_vec(2,1,:), [], 1);
            A22 = reshape(A_vec(2,2,:), [], 1);

            % build global row/col/value arrays (4 entries per element)
            %   i.e. in effect, for each element e with nodes [i j], we add:
            %   (i,i) += A11(e)
            %   (i,j) += A12(e)
            %   (j,i) += A21(e)
            %   (j,j) += A22(e)
            I = [ elems(:,1); elems(:,1); elems(:,2); elems(:,2) ];  % 4*nElems x 1
            J = [ elems(:,1); elems(:,2); elems(:,1); elems(:,2) ];  % 4*nElems x 1
            V = [ A11;        A12;        A21;        A22        ];  % 4*nElems x 1

            % final assembly
            A = sparse(I, J, V, nNodes, nNodes);

        end

        function b = computeVolumeForces(self, f1, f2)
        % COMPUTEVOLUMEFORCES Vectorized computation of volume load vector.
        % 
        % Usage:
        %   b = computeVolumeForces(self)
        %   b = computeVolumeForces(self, f1)
        %   b = computeVolumeForces(self, f1, f2)   
        % Inputs:
        %   f1, f2 : (optional) function handles for volume force factors
        %            after checkVariables, both are used as f(x,t,u,v)
        %            If not provided, f1 = self.f, f2 ≡ 1.
        % Output:
        %   b : nNodes x 1 sparse load vector

            % input handling
            if nargin == 1
                f1 = self.f;
                f2 = @(x,t) 1;
            elseif nargin == 2
                f2 = @(x,t) 1;
            end

            % unpack mesh avariables
            coords = self.domain.mesh.nodes(:);
            nNodes = numel(coords);
            elems  = self.domain.mesh.elements;
            nElems = size(elems, 1);
            mids   = self.domain.mesh.midpoints;
            areas  = self.domain.mesh.areas;

            % check function signatures, i.e. f = f(x,t,U,V)
            [f1, t, U, V] = self.checkVariables(f1);
            f2 = self.checkVariables(f2);

            % Interpolate U to midpoints 
            U_mid = mean(U(elems), 2);   

            % compute volume forces
            f1_vals = f1(mids,t,U_mid,V); 
            f2_vals = f2(mids,t,U_mid,V); 
            forces  = f1_vals .* f2_vals / 2; 
            volumeForces = areas .* forces;

            % each element e contributes volumeForces(e) to each of its 2 nodes
            b = accumarray(...
                    elems(:), ...
                    repmat(volumeForces, 2, 1), ...
                    [nNodes, 1], ...
                    @sum, ...
                    0, ...
                    true);

            % make sure it's sparse (accumarray with 'true' already does this)
            b = sparse(b);
        end
        

        % BOUNDARY CONDITION FUNCTIONS
        %{
		function self = assembleBCs(self)

			%self = self.computePeriodicBCs;
			self.vectors.U_D   = self.computeDirichletBCs;
			%self.vectors.b_neu = self.computeNeumannBCs;
			%[self.tensors.M_rob,self.vectors.b_rob] = self.computeRobinBCs;

		end
        %}

        function U_D = computeDirichletBCs(self)

            % unpack variables
            dom    = self.domain;
            nNodes = dom.mesh.nNodes;
            coords = dom.mesh.nodes(:);   % ensure column vector

            % initialize storage
            U_D = sparse(nNodes, 1);

            % loop over boundary edges
            for i = 1:dom.boundary.nEdges

                edge = dom.boundary.edges(i);

                if edge.boundaryType == 'D'

                    % nodes on i-th boundary edge
                    bNodes_i = edge.nodes;

                    % grab current BC handle
                    bc = edge.boundaryCondition;

                    % check if boundary condition is time-varying
                    if Coefficients.isTimeVarying(bc) == 0
                        % wrap into (x,t) form, with t ignored
                        bc = @(x,t) bc(x);
                        t  = 0;
                    else
                        t  = self.t;
                    end

                    % compute Dirichlet values at boundary nodes
                    U_D(bNodes_i) = bc(coords(bNodes_i), t);

                end
            end
        end

        function b_neu = computeNeumannBCs(self)

            % unpack variables
            dom    = self.domain;
            nNodes = dom.mesh.nNodes;
            coords = dom.mesh.nodes(:);   % ensure column

            % initialize storage
            b_neu = sparse(nNodes, 1);

            % loop over boundary edges
            for i = 1:dom.boundary.nEdges

                edge = dom.boundary.edges(i);

                % only process Neumann edges
                if edge.boundaryType == 'N'

                    % nodes on this Neumann boundary edge
                    bNodes_i = edge.nodes;          % indices into global node list
                    bCond    = edge.boundaryCondition;

                    % normalize boundary condition into f(x,t,u,v) form
                    [bCond, t, U, V] = self.checkVariables(bCond);

                    % coordinates and solution values at these boundary nodes
                    x_b   = coords(bNodes_i);
                    u_b   = U(bNodes_i);

                    % evaluate Neumann flux at boundary nodes
                    q = bCond(x_b, t, u_b, V);      % could be scalar or vector

                    if numel(q) == 1
                        q = q * ones(numel(bNodes_i), 1);  % broadcast scalar
                    end

                    % accumulate into global RHS
                    b_neu(bNodes_i) = b_neu(bNodes_i) + q;
                end
            end
        end

        function [E, b] = computeRobinBCs(self)

            % unpack variables
            dom    = self.domain;
            nNodes = dom.mesh.nNodes;
            coords = dom.mesh.nodes(:);   % ensure column

            % initialize storage
            b = sparse(nNodes, 1);
            E = sparse(nNodes, nNodes);

            % loop over boundary edges
            for i = 1:dom.boundary.nEdges

                edge = dom.boundary.edges(i);

                % only process Robin edges
                if edge.boundaryType == 'R'

                    % boundaryCondition is assumed to be {alpha, u_BC}
                    bCond = edge.boundaryCondition;

                    alpha = bCond{1};
                    u_BC  = bCond{2};

                    % normalize both into f(x,t,u,v) form
                    [alpha, t, U, V] = self.checkVariables(alpha);
                    [u_BC, t, U, V]  = self.checkVariables(u_BC);

                    % nodes on this Robin boundary edge
                    bNodes_i = edge.nodes;     % indices of boundary nodes
                    x_b      = coords(bNodes_i);
                    u_b      = U(bNodes_i);

                    % evaluate alpha and u_BC at boundary nodes
                    alpha_vals = alpha(x_b, t, u_b, V);   % may be scalar or vector
                    uBC_vals   = u_BC(x_b, t, u_b, V);

                    % broadcast scalar values if needed
                    if numel(alpha_vals) == 1
                        alpha_vals = alpha_vals * ones(numel(bNodes_i), 1);
                    end
                    if numel(uBC_vals) == 1
                        uBC_vals = uBC_vals * ones(numel(bNodes_i), 1);
                    end

                    % RHS contribution: alpha * u_BC
                    b(bNodes_i) = b(bNodes_i) + alpha_vals .* uBC_vals;

                    % Matrix contribution: alpha * u
                    % (purely diagonal, one per Robin node)
                    E = E + sparse(bNodes_i, bNodes_i, alpha_vals, nNodes, nNodes);
                end
            end
        end

        % UTILITY FUNCTIONS
        function [f,t,U,V] = checkVariables(self,f)

            % store function type data
            isTimeVarying = Coefficients1d.isTimeVarying(f);    
            isNonlinear   = Coefficients1d.isNonlinear(f);
            isCoupled     = Coefficients1d.isCoupled(f);
            
            % expand signature to f = f(x,t,u,v)
            f = Coefficients1d.checkVariables(f);

            % prepare t, U, V
            if isTimeVarying
                t = self.t;
            else
                t = 0;
            end

            if isNonlinear 
                U = self.U;
            else
                U = zeros(self.domain.mesh.nNodes,1);
            end

            if isCoupled
                V = self.V;
            else
                V = zeros(self.domain.mesh.nNodes,1);
            end

        end
    
		function self = cleanup(self)

            % clear out tensors and vectors
			self.tensors = [];
			self.vectors = [];

			%{
			% copy solution to periodic replica nodes
			P = self.domain.boundary.P_nodes;
			self.solution(P.replica.edge,:) = self.solution(P.free.edge,:);
			self.solution(P.replica.corner,:) = ...
							repmat(self.solution(P.free.corner,:),3,1);					
			%self.solution(P.replica.corner,:) = self.solution(P.free.corner,:);
			%}

            % ensure solution is full
            self.solution = full(self.solution);

		end


        % PLOTTERS
        function plot(self,timestep)

            % input handling
            if nargin < 2 && ~isempty(self.domain.time)
                timestep = self.domain.time.N_T;

            elseif nargin < 2 && isempty(self.domain.time)
                timestep = 1;
            end

            % use Domain1d's plot function
            sol = self.solution(:,timestep);
            self.domain.plot_function(sol);

            %{
            % extract time and solution at requested timestep
            t = self.domain.time.timeSteps(timestep);
            u = self.solution(:, timestep);

            % plot
            figure;
            plot(self.domain.mesh.nodes, u, '-o', 'LineWidth', 2);
            title(['Solution at t = ', num2str(t)]);
            xlabel('x');
            ylabel('u(x,t)');
            grid on;
            %}


        end


    end
end


