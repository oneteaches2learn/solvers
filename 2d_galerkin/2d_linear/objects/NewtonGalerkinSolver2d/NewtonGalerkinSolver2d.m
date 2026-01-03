classdef NewtonGalerkinSolver2d

    properties
    end

    methods
        function self = NewtonGalerkinSolver2d()
            ...
        end

        function E = computeNonlinearNeumannContribution(self)

			switch self.options.assemblyQuadrature
				case 'legacy'
					E = self.computeNonlinearNeumannContribution_legacy();
				case 'centroid'
					E = self.computeNonlinearNeumannContribution_centroid();
				case 'gaussP1'
					E = self.computeNonlinearNeumannContribution_gaussP1();
				otherwise
					error('Unknown assemblyQuadrature option.');
			end

		end



        function E = computeNonlinearNeumannContribution_legacy(self)

			% unpack variables
			dom    = self.domain;
			nNodes = self.domain.mesh.nNodes;
			coords = self.domain.mesh.nodes;

			% initialize storage
			E = sparse(nNodes,nNodes);

			% compute boundary conditions
			for i = 1:self.domain.boundary.nEdges
				
				% compute Dirichlet condition
				if dom.boundary.edges(i).boundaryType == 'N'
					
					bCond = dom.boundary.edges(i).boundaryCondition_ddu;

                    % check if alpha is time-varying
					[bCond,t,U,V] = self.checkVariables(bCond);
                    
					% store nodes on i-th edge of domain
					bNodes_i = dom.boundary.edges(i).nodes;

					% loop over segments of i-th edge
					for j = 1:length(bNodes_i)-1

						% get edge data
						edge = [bNodes_i(j) bNodes_i(j+1)];
						edgeMidPt = sum(coords(edge,:)/2);
						edgeLength = norm(coords(edge(1),:) - coords(edge(2),:));

						% compute E matrix
						E(edge(1),edge(1)) = E(edge(1),edge(1)) + ...
							1/2 * edgeLength * bCond(coords(edge(1),1),coords(edge(1),2),t,U(edge(1)),V);
						E(edge(2),edge(2)) = E(edge(2),edge(2)) + ...
							1/2 * edgeLength * bCond(coords(edge(2),1),coords(edge(2),2),t,U(edge(2)),V);
					end
				end
			end
		end
		% =====================================================================

	%{
% this is an older version to be deleted possibly
function E = computeNonlinearNeumannContribution_centroid(self)
% Vectorized centroid (midpoint) assembly for nonlinear Neumann DJ term:
%   E_ij += ∫_{Gamma_N} g'(u_h) phi_i phi_j ds
% using midpoint quadrature on each boundary segment.
%
% Per-edge safe: allows different Neumann edges to have different bCond handles.

    dom    = self.domain;
    nNodes = dom.mesh.nNodes;
    coords = dom.mesh.nodes;

    % We'll accumulate COO-format triplets and build sparse once.
    Iall = [];
    Jall = [];
    Sall = [];

    for i = 1:dom.boundary.nEdges

        if dom.boundary.edges(i).boundaryType ~= 'N'
            continue;
        end

        % derivative handle on this edge
        bCond = dom.boundary.edges(i).boundaryCondition_ddu;
        [bCond,t,U,V] = self.checkVariables(bCond);

        % nodes on this boundary polyline
        bn = dom.boundary.edges(i).nodes(:);
        if numel(bn) < 2
            continue;
        end

        % segments on this edge (vectorized over segments)
        e1 = bn(1:end-1);
        e2 = bn(2:end);

        x1 = coords(e1,:);
        x2 = coords(e2,:);
        xm = 0.5*(x1 + x2);

        dx = x2 - x1;
        L  = sqrt(sum(dx.^2,2));

        Um = 0.5*(U(e1) + U(e2));

        % evaluate coefficient c = g'(u_h) at midpoints
        % (try vectorized; fallback to scalar calls)
        try
            c = bCond(xm(:,1), xm(:,2), t, Um, V);
        catch
            nSeg = size(xm,1);
            c = arrayfun(@(k) bCond(xm(k,1), xm(k,2), t, Um(k), V), (1:nSeg).');
        end

        a = (L .* c) / 6;   % midpoint weight

        % Add the 2x2 local matrix (L*c/6)*[2 1; 1 2] for every segment
        I = [e1; e1; e2; e2];
        J = [e1; e2; e1; e2];
        S = [2*a; 1*a; 1*a; 2*a];

        Iall = [Iall; I];
        Jall = [Jall; J];
        Sall = [Sall; S];
    end

    if isempty(Iall)
        E = sparse(nNodes,nNodes);
        return;
    end

    E = sparse(Iall, Jall, Sall, nNodes, nNodes);
end
%}


function E = computeNonlinearNeumannContribution_centroid(self)
% Clean, per-edge safe, vectorized CENTROID (midpoint) assembly for nonlinear Neumann DJ term:
%   E_ij += ∫_{Gamma_N} g'(u_h) phi_i phi_j ds
% using midpoint quadrature on each boundary segment.
%
% Assumes boundaryCondition_ddu is VECTORIZED in (x,y,u): accepts vectors and returns one value per point.

    dom    = self.domain;
    nNodes = dom.mesh.nNodes;
    coords = dom.mesh.nodes;

    % Accumulate COO triplets and build sparse once.
    Iall = [];
    Jall = [];
    Sall = [];

    for i = 1:dom.boundary.nEdges

        if dom.boundary.edges(i).boundaryType ~= 'N'
            continue;
        end

        % edge-specific derivative handle
        bCond = dom.boundary.edges(i).boundaryCondition_ddu;
        [bCond,t,U,V] = self.checkVariables(bCond);

        bn = dom.boundary.edges(i).nodes(:);
        if numel(bn) < 2
            continue;
        end

        % segments (vectorized over segments on this boundary edge)
        e1 = bn(1:end-1);
        e2 = bn(2:end);

        x1 = coords(e1,:);
        x2 = coords(e2,:);
        xm = 0.5*(x1 + x2);

        dx = x2 - x1;
        L  = sqrt(sum(dx.^2,2));

        Um = 0.5*(U(e1) + U(e2));

        % vectorized coefficient eval at midpoints
        c = bCond(xm(:,1), xm(:,2), t, Um, V);
        c = c(:);

        % optional sanity check (uncomment if you like)
        % assert(numel(c)==numel(Um), 'boundaryCondition_ddu must return one value per midpoint.');

        a = (L .* c) / 6;   % (nSeg x 1)

        % add (L*c/6)*[2 1; 1 2] per segment
        I = [e1; e1; e2; e2];
        J = [e1; e2; e1; e2];
        S = [2*a; 1*a; 1*a; 2*a];

        Iall = [Iall; I];
        Jall = [Jall; J];
        Sall = [Sall; S];
    end

    if isempty(Iall)
        E = sparse(nNodes,nNodes);
        return;
    end

    E = sparse(Iall, Jall, Sall, nNodes, nNodes);
end


		%{
		function E = computeNonlinearNeumannContribution_centroid(self)

			% unpack variables
			dom    = self.domain;
			nNodes = dom.mesh.nNodes;
			coords = dom.mesh.nodes;

			% initialize storage
			E = sparse(nNodes,nNodes);

			% loop over boundary edges
			for i = 1:dom.boundary.nEdges

				if dom.boundary.edges(i).boundaryType ~= 'N'
					continue;
				end

				% this is dg/du (or whatever you store as the Newton derivative)
				bCond = dom.boundary.edges(i).boundaryCondition_ddu;

				% check variable conventions (time/state dependence)
				[bCond,t,U,V] = self.checkVariables(bCond);

				% nodes on this boundary polyline
				bNodes_i = dom.boundary.edges(i).nodes;

				% loop over segments
				for j = 1:length(bNodes_i)-1

					edge = [bNodes_i(j) bNodes_i(j+1)];

					x1 = coords(edge(1),:);
					x2 = coords(edge(2),:);
					L  = norm(x1 - x2);

					% midpoint (centroid on segment)
					xm = 0.5*(x1 + x2);

					% P1 midpoint value
					Um = 0.5*(U(edge(1)) + U(edge(2)));

					% coefficient at midpoint
					c_m = bCond(xm(1), xm(2), t, Um, V);

					% local edge "mass" matrix: (L*c_m/6)*[2 1; 1 2]
					E(edge,edge) = E(edge,edge) + (L * c_m / 6) * [2 1; 1 2];

				end
			end
		end
		%}

		%{
		function E = computeNonlinearNeumannContribution_gaussP1(self)

			% unpack variables
			dom    = self.domain;
			nNodes = dom.mesh.nNodes;
			coords = dom.mesh.nodes;

			% initialize storage
			E = sparse(nNodes,nNodes);

			% 2-point Gauss rule on [0,1]
			s = [0.5 - 1/(2*sqrt(3)), 0.5 + 1/(2*sqrt(3))];
			w = [0.5, 0.5];

			% loop over boundary edges
			for i = 1:dom.boundary.nEdges

				if dom.boundary.edges(i).boundaryType ~= 'N'
					continue;
				end

				% this is dg/du (Newton derivative contribution)
				bCond = dom.boundary.edges(i).boundaryCondition_ddu;

				% check variable conventions (time/state dependence)
				[bCond,t,U,V] = self.checkVariables(bCond);

				% nodes on this boundary polyline
				bNodes_i = dom.boundary.edges(i).nodes;

				% loop over segments
				for j = 1:length(bNodes_i)-1

					edge = [bNodes_i(j) bNodes_i(j+1)];

					x1 = coords(edge(1),:);
					x2 = coords(edge(2),:);
					U1 = U(edge(1));
					U2 = U(edge(2));

					L  = norm(x1 - x2);

					E_loc = zeros(2,2);

					for q = 1:2
						sq  = s(q);
						phi = [1 - sq; sq];          % [phi1; phi2] on the segment
						xq  = (1 - sq)*x1 + sq*x2;   % point on segment
						Uq  = phi(1)*U1 + phi(2)*U2; % P1 exact on the segment

						cq  = bCond(xq(1), xq(2), t, Uq, V);

						E_loc = E_loc + (L * w(q) * cq) * (phi*phi.');
					end

					E(edge,edge) = E(edge,edge) + E_loc;

				end
			end
		end
		% =====================================================================
		%}

function E = computeNonlinearNeumannContribution_gaussP1(self)
% Vectorized 2-point Gauss assembly for nonlinear Neumann DJ term:
%   E_ij += ∫_{Gamma_N} g'(u_h) phi_i phi_j ds
% using 2-point Gauss quadrature on each boundary segment.
%
% Per-edge safe: allows different Neumann edges to have different bCond handles.

    dom    = self.domain;
    nNodes = dom.mesh.nNodes;
    coords = dom.mesh.nodes;

    % 2-point Gauss rule on [0,1]
    s = [0.5 - 1/(2*sqrt(3)), 0.5 + 1/(2*sqrt(3))];
    w = [0.5, 0.5];

    % P1 basis values at Gauss points on [0,1]
    phi1 = 1 - s;   % row vector length 2
    phi2 = s;

    % We'll accumulate COO-format triplets and build sparse once.
    Iall = [];
    Jall = [];
    Sall = [];

    for i = 1:dom.boundary.nEdges

        if dom.boundary.edges(i).boundaryType ~= 'N'
            continue;
        end

        % derivative handle on this edge
        bCond = dom.boundary.edges(i).boundaryCondition_ddu;
        [bCond,t,U,V] = self.checkVariables(bCond);

        % nodes on this boundary polyline
        bn = dom.boundary.edges(i).nodes(:);
        if numel(bn) < 2
            continue;
        end

        % segments on this edge
        e1 = bn(1:end-1);
        e2 = bn(2:end);

        x1 = coords(e1,:);     % (nSeg x 2)
        x2 = coords(e2,:);
        dx = x2 - x1;
        L  = sqrt(sum(dx.^2,2));  % (nSeg x 1)

        U1 = U(e1);
        U2 = U(e2);

        % --- Gauss point 1 (vectorized over segments)
        xq1 = phi1(1).*x1 + phi2(1).*x2;     % (nSeg x 2)
        Uq1 = phi1(1).*U1 + phi2(1).*U2;     % (nSeg x 1)

        % --- Gauss point 2
        xq2 = phi1(2).*x1 + phi2(2).*x2;
        Uq2 = phi1(2).*U1 + phi2(2).*U2;

        % Evaluate coefficient c = g'(u_h) at Gauss points
        % Try vectorized first; fallback otherwise
        try
            c1 = bCond(xq1(:,1), xq1(:,2), t, Uq1, V);
            c2 = bCond(xq2(:,1), xq2(:,2), t, Uq2, V);
        catch
            nSeg = size(x1,1);
            c1 = arrayfun(@(k) bCond(xq1(k,1), xq1(k,2), t, Uq1(k), V), (1:nSeg).');
            c2 = arrayfun(@(k) bCond(xq2(k,1), xq2(k,2), t, Uq2(k), V), (1:nSeg).');
        end

        % Local 2x2 matrix entries per segment:
        % E_loc += sum_q (L*w_q*c_q) * [phi1_q;phi2_q]*[phi1_q phi2_q]
        %
        % Compute per-segment scalars for each entry:
        a11 = L .* ( w(1)*c1*phi1(1)*phi1(1) + w(2)*c2*phi1(2)*phi1(2) );
        a12 = L .* ( w(1)*c1*phi1(1)*phi2(1) + w(2)*c2*phi1(2)*phi2(2) );
        a22 = L .* ( w(1)*c1*phi2(1)*phi2(1) + w(2)*c2*phi2(2)*phi2(2) );

        % symmetry: a21 = a12
        I = [e1; e1; e2; e2];
        J = [e1; e2; e1; e2];
        S = [a11; a12; a12; a22];

        Iall = [Iall; I];
        Jall = [Jall; J];
        Sall = [Sall; S];
    end

    if isempty(Iall)
        E = sparse(nNodes,nNodes);
        return;
    end

    E = sparse(Iall, Jall, Sall, nNodes, nNodes);
end


		function E = computeNonlinearRobinContribution(self)

			% unpack variables
			dom    = self.domain;
			nNodes = self.domain.mesh.nNodes;
			coords = self.domain.mesh.nodes;

			% initialize storage
			E = sparse(nNodes,nNodes);

			% compute boundary conditions
			for i = 1:self.domain.boundary.nEdges
				
				% compute Dirichlet condition
				if dom.boundary.edges(i).boundaryType == 'R'
					
					alpha = dom.boundary.edges(i).boundaryCondition{1};
					bCond = dom.boundary.edges(i).boundaryCondition_ddu;

                    % check if alpha is time-varying
					[bCond,t,U,V] = self.checkVariables(bCond);
					[alpha,t,U,V] = self.checkVariables(alpha);
                    
					% store nodes on i-th edge of domain
					bNodes_i = dom.boundary.edges(i).nodes;

					% loop over segments of i-th edge
					for j = 1:length(bNodes_i)-1

						% get edge data
						edge = [bNodes_i(j) bNodes_i(j+1)];
						edgeMidPt = sum(coords(edge,:)/2);
						edgeLength = norm(coords(edge(1),:) - coords(edge(2),:));

						% compute E matrix
						E(edge(1),edge(1)) = E(edge(1),edge(1)) + ...
							1/2 * edgeLength * bCond(coords(edge(1),1),coords(edge(1),2),t,sum(U(edge))/2,V);
						E(edge(2),edge(2)) = E(edge(2),edge(2)) + ...
							1/2 * edgeLength * bCond(coords(edge(2),1),coords(edge(2),2),t,sum(U(edge))/2,V);

					end
				end
			end
		end

    end
end

    