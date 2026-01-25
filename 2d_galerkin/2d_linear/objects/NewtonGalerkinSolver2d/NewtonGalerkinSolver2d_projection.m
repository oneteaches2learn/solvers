classdef NewtonGalerkinSolver2d_projection < NewtonGalerkinSolver2d_elliptic
% NewtonGalerkinSolver2d_projection
%
% Computes W = P_h(u;v) for fixed v by solving:
%   (k∇W,∇z) + lambda(W,z) - (alpha(W,v),z)_∂Ω
%     = (k∇u,∇z) + lambda(u,z) - (alpha(u,v),z)_∂Ω
%
% Implementation notes:
%   - LHS nonlinear boundary term is implemented via boundaryType 'N' using
%     boundaryCondition = @(x,y,t,U,V) -alpha(U,V)
%   - RHS uses u_fixed, v_fixed to build:
%       b_vol_proj = (A + M_lambda) * u_vec
%       b_alpha_u  = boundary vector for alpha(u_fixed,v_fixed)

    properties
        u_vec          % nodal interpolation of u_fixed at mesh nodes
        b_alpha_u      % RHS boundary vector for (alpha(u_fixed,v_fixed), z)_∂Ω
        V_fixed        % fixed parameter v
        U_fixed        % nodal interpolation of u_fixed
    end

    properties (Hidden)
        V
        V_prev         % V_prev is used by coupled solver for computations involving v
    end

    methods
        function self = NewtonGalerkinSolver2d_projection(dom, auxfun)

            % call superclass constructor (stores dom, cofs, f)
            self@NewtonGalerkinSolver2d_elliptic(dom, auxfun);

            % store projection parameters
            self.coefficients.u_fixed = auxfun.params.u_fixed;
            self.coefficients.v_fixed = auxfun.params.v_fixed;
            self.coefficients.lambda  = auxfun.params.lambda;

			% evaluate projected function on mesh nodes
			nodes = self.domain.mesh.nodes;
			U_fixed = self.coefficients.u_fixed(nodes(:,1),nodes(:,2));
			V_fixed = self.coefficients.v_fixed;

			% set initial guess to (U_fixed,V_fixed)
			self.U = U_fixed;
            self.U_fixed = U_fixed;
			self.V = V_fixed;
            self.V_fixed = V_fixed;
            self.V_prev = V_fixed;
			%self.U = ones(size(self.U)) * mean(U_fixed);  % better initial guess
			%self.U = ones(size(self.U));

        end

        function self = assembleVectors(self)
            % We override assembleVectors because we do NOT want
            % b_vol = computeVolumeForces(self.f). Instead we want
            % b_vol = (A + M_lambda)*u_vec minus boundary alpha(u,v) term.

            % assemble vector
			self.vectors.M_r = self.computeVolumeForces(self.coefficients.r);

            % Ensure tensors are already assembled when this is called.
            % (Your solve() calls assembleTensors before assembleVectors.)
            A   = self.tensors.A;
            M_dr = self.tensors.M_dr;  
            lambda = self.coefficients.lambda;

            % Volume part of RHS: (k∇u,∇z)+lambda(u,z)  -> (A+M_r)*u_vec
            %self.vectors.b_vol = (A + M_dr) * self.U_fixed;
            self.vectors.b_vol = (A + lambda * M_dr) * self.U_fixed;

            % Boundary RHS part: -(alpha(u_fixed,v_fixed), z)_∂Ω.
            % We will assemble b_alpha_u := ∫ alpha(u_fixed,v_fixed) z  (same structure as Neumann load)
            % and then subtract it in finalAssembly.
            self.vectors.b_alpha_u = self.computeAlphaUVector();
        end

        function [DJ, J] = finalAssembly(self, U_tilde)
            % Matches your NewtonGalerkinSolver2d_poisson pattern, but:
            %   - includes M_r in the linear operator S (since lambda is linear)
            %   - subtracts RHS boundary term b_alpha_u

            %{
            tensors = self.tensors;
            vectors = self.vectors;

            % Linear operator
            S = tensors.A + tensors.M_r + tensors.M_rob;

            % Load vector (RHS)
            % NOTE: b_neu is computed from current iterate U via computeNeumannBCs,
            % and is SUBTRACTED here, so it will appear with a PLUS sign in the residual J,
            % exactly as in your existing Newton formulation.
            %
            % We also subtract b_alpha_u to enforce the RHS term -(alpha(u,v),z)_∂Ω.
            b = (vectors.b_vol ...
                 - self.b_alpha_u ...
                 - vectors.b_neu ...
                 + vectors.b_rob) ...
                 - S * vectors.U_D;

            % Residual
            J = S * U_tilde - b;

            % Jacobian
            DJ = tensors.A + tensors.M_r + tensors.M_rob ...
                 + tensors.M_dr + tensors.M_dneu - tensors.M_drob;
            %}

			% store variables
            tensors = self.tensors;
            vectors = self.vectors;

			% assemble Linear Tensor
			S = tensors.A;

			% assemble nonlinear contributions to J
			b_nonlinear = self.vectors.M_r;

			% assemble Load Vector
			b = (vectors.b_vol + vectors.b_neu - vectors.b_alpha_u);

			% assemble J
	        J = S * U_tilde + b_nonlinear - b;
	        J = (S + tensors.M_dr) * U_tilde - b;
			%J = (tensors.A + tensors.M_dr) * U_tilde - b;
			%J = (tensors.A + tensors.M_dr - tensors.M_dneu) * U_tilde - vectors.b_vol;
                
			% assemble DJ
			DJ = tensors.A + tensors.M_dr - tensors.M_dneu;

        end
    end

    methods (Access = private)
        function b_alpha = computeAlphaUVector(self)
            % Assemble boundary vector b_alpha with integrand alpha(u_fixed(x,y), v_fixed)
            % using the same midpoint-per-segment structure as computeNeumannBCs.
            %
            % This produces b_alpha ≈ ∫_{∂Ω_N} alpha(u_fixed,v_fixed) * phi_i ds.

            % unpack data
            dom    = self.domain;
            nNodes = dom.mesh.nNodes;
            coords = dom.mesh.nodes;
            v_fixed = self.V_fixed;
            u_fixed = self.coefficients.u_fixed;
            t = 0; % placeholder, since problem is time-independent

            % assemble b_alpha vector
            b_alpha = sparse(nNodes,1);

            for i = 1:dom.boundary.nEdges
                if dom.boundary.edges(i).boundaryType ~= 'N'
                    continue;
                end

                % store computational data
                bn = dom.boundary.edges(i).nodes;
                bCond = dom.boundary.edges(i).boundaryCondition;
				[bCond,t,U,V] = self.checkVariables(bCond);

                for j = 1:length(bn)-1
                    edge = [bn(j) bn(j+1)];
                    xm = sum(coords(edge,:),1)/2;
                    L  = norm(coords(edge(1),:) - coords(edge(2),:));

                    u_mid = u_fixed(xm(1), xm(2));
                    val   = bCond(xm(1), xm(2), t, u_mid, v_fixed);

                    % consistent with computeNeumannBCs midpoint rule:
                    b_alpha(edge) = b_alpha(edge) + (L * val / 2);
                end
            end

        end
    end
end
