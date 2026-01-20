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
    end

    methods
        function self = NewtonGalerkinSolver2d_projection(dom, auxfun)

            % call superclass constructor (stores dom, cofs, f)
            self@NewtonGalerkinSolver2d_elliptic(dom, auxfun);

            % solve immediately (matches your other subclasses)
            % self = self.solve;

        end

        function self = assembleVectors(self)
            % We override assembleVectors because we do NOT want
            % b_vol = computeVolumeForces(self.f). Instead we want
            % b_vol = (A + M_lambda)*u_vec minus boundary alpha(u,v) term.

            % Ensure tensors are already assembled when this is called.
            % (Your solve() calls assembleTensors before assembleVectors.)
            A   = self.tensors.A;
            M_r = self.tensors.M_r;   % this should be lambda mass matrix if cofs.r = lambda

            % Evaluate u_fixed at nodes to get nodal interpolation u_h.
            coords = self.domain.mesh.nodes;
            if ~isfield(self, 'f') && ~isfield(self, 'coefficients')
                error('Projection solver: unexpected missing fields. Check construction.');
            end
            if ~isfield(self, 'domain')
                error('Projection solver: domain not set.');
            end

            % auxfun.u_fixed is not stored directly; it lives in self.f only if you set it.
            % Easiest: store u_fixed in auxfun.f and pass it in as auxfun.f = @(x,y,...) ...
            % But since you asked for auxfun.u_fixed, we will fetch it from coefficients if you store it there.
            %
            % Recommended: add u_fixed, v_fixed, alpha, dalpha_du to auxfun itself
            % AND ALSO store u_fixed on self as self.f? (see usage snippet below).
            %
            % Here we assume you placed u_fixed in self.domain.userdata or self.options.
            % If you follow the usage snippet below, we will store it in self.options.projection.
            if ~isfield(self.options, 'projection') || ~isfield(self.options.projection, 'u_fixed')
                error(['Projection solver needs self.options.projection.u_fixed. ', ...
                       'See usage snippet below.']);
            end

            u_fixed = self.options.projection.u_fixed;

            x = coords(:,1); y = coords(:,2);
            self.u_vec = u_fixed(x,y);

            % Volume part of RHS: (k∇u,∇z)+lambda(u,z)  -> (A+M_r)*u_vec
            self.vectors.b_vol = (A + M_r) * self.u_vec;

            % Boundary RHS part: -(alpha(u_fixed,v_fixed), z)_∂Ω.
            % We will assemble b_alpha_u := ∫ alpha(u_fixed,v_fixed) z  (same structure as Neumann load)
            % and then subtract it in finalAssembly.
            self.b_alpha_u = self.computeAlphaUVector();
        end

        function [DJ, J] = finalAssembly(self, U_tilde)
            % Matches your NewtonGalerkinSolver2d_poisson pattern, but:
            %   - includes M_r in the linear operator S (since lambda is linear)
            %   - subtracts RHS boundary term b_alpha_u

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
        end
    end

    methods (Access = private)
        function b_alpha = computeAlphaUVector(self)
            % Assemble boundary vector b_alpha with integrand alpha(u_fixed(x,y), v_fixed)
            % using the same midpoint-per-segment structure as computeNeumannBCs.
            %
            % This produces b_alpha ≈ ∫_{∂Ω_N} alpha(u_fixed,v_fixed) * phi_i ds.

            dom    = self.domain;
            nNodes = dom.mesh.nNodes;
            coords = dom.mesh.nodes;

            if ~isfield(self.options, 'projection') || ~isfield(self.options.projection, 'alpha')
                error(['Projection solver needs self.options.projection.alpha ', ...
                       '(handle alpha(u,v)). See usage snippet below.']);
            end
            if ~isfield(self.options.projection, 'v_fixed')
                error('Projection solver needs self.options.projection.v_fixed.');
            end
            if ~isfield(self.options.projection, 'u_fixed')
                error('Projection solver needs self.options.projection.u_fixed.');
            end

            alpha  = self.options.projection.alpha;
            v_fixed = self.options.projection.v_fixed;
            u_fixed = self.options.projection.u_fixed;

            b_alpha = sparse(nNodes,1);

            for i = 1:dom.boundary.nEdges
                if dom.boundary.edges(i).boundaryType ~= 'N'
                    continue;
                end

                bn = dom.boundary.edges(i).nodes;

                for j = 1:length(bn)-1
                    edge = [bn(j) bn(j+1)];
                    xm = sum(coords(edge,:),1)/2;
                    L  = norm(coords(edge(1),:) - coords(edge(2),:));

                    u_mid = u_fixed(xm(1), xm(2));
                    val   = alpha(u_mid, v_fixed);

                    % consistent with computeNeumannBCs midpoint rule:
                    b_alpha(edge) = b_alpha(edge) + (L * val / 2);
                end
            end
        end
    end
end
