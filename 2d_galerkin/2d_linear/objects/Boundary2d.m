classdef Boundary2d

	properties
		nEdges
		edges
		D_nodes
		N_nodes
		R_nodes
		P_nodes
		freeNodes
	end

	methods
		% CONSTUCTOR
		function self = Boundary2d(dl,bcTypes,bcConds)

			% create edge geometry
			if nargin > 0
				self.nEdges = self.set_nEdges(dl);
				self.edges = self.setEdgeGeometry(dl);
			end

			% set boundary condition types
			if nargin > 1
				self = self.setBCTypes(bcTypes);
			end

			% set boundary conditions
			if nargin > 2
				self = self.setBCConditions(bcConds);
			end

		end


		% SETTERS
		function nEdges = set_nEdges(self,dl)

			nEdges = size(dl,2);

		end

		function edges = setEdgeGeometry(self,dl)

			edges(self.nEdges) = BoundaryEdge2d_test;
			for j = 1:self.nEdges;

				% create new edge
				vert = reshape(dl(2:5,j),2,2);
				edges(j).vertex1 = vert(1,:);
				edges(j).vertex2 = vert(2,:);
				edges(j).ID = j;

			end

		end

		function self = setBCTypes(self,bcTypes)

			% convert cell array to char array if necessary
			if isa(bcTypes,'cell')
				for i = length(bcTypes)
					temp(i) = bcTypes{i};
				end
				bcTypes = temp;
			end

			% if a bcType is given for each edge
			if length(bcTypes) == self.nEdges
				for i = 1:self.nEdges
					self.edges(i).boundaryType = bcTypes(i);
				end

			% if one bcType is prescribed for all internal edges
			else
				for i = 1:4
					self.edges(i).boundaryType = bcTypes(i);
				end
				for i = 5:self.nEdges
					self.edges(i).boundaryType = bcTypes(5);
					bcTypes(i) = bcTypes(5);
				end
			end

		end

		function self = setBCConditions(self,bcConds)

			for i = 1:self.nEdges
				self.edges(i).boundaryCondition = bcConds{i};
			end

		end

		function self = setEdgeNodes(self,mesh)

			for i = 1:self.nEdges
				edgeID = self.edges(i).ID;
				self.edges(i).nodes = mesh.Mesh.findNodes('region','Edge',edgeID);
				self.edges(i).nNodes = length(self.edges(i).nodes);
			end

		end

		function self = setFreeNodes(self,mesh)

			boundNodes = [self.D_nodes, self.P_nodes.replica.edge, ...
												self.P_nodes.replica.corner];
			self.freeNodes = setdiff(1:mesh.nNodes,boundNodes);

		end

		function self = setBoundaryNodeLists(self,mesh)

			% initialize boundary node lists
			D = [];
			N = [];
			R = [];
			P = [];

			% compile boundary node lists
			for i = 1:self.nEdges
				if strcmp(self.edges(i).boundaryType,'D')
					D = [D self.edges(i).nodes];
				elseif strcmp(self.edges(i).boundaryType,'N')
					N = [N self.edges(i).nodes];
				elseif strcmp(self.edges(i).boundaryType,'R')
					R = [R self.edges(i).nodes];
				end
			end
			P = self.setPeriodicNodes;

			% sort D nodes and remove duplicates
			D = unique(D);

			% remove D nodes from R and N node lists, sort and remove duplicates
			N = setdiff(N,D);
			R = setdiff(R,D);

			% store boundary node lists
			self.D_nodes = D;
			self.N_nodes = N;
			self.R_nodes = R;
			self.P_nodes = P;

		end

		function P = setPeriodicNodes(self)

			% build periodic node dictionary
			P.free.edge = [];
			P.free.corner = [];
			P.replica.edge = [];
			P.replica.corner = [];
			
			% if North-South periodicity
			if (self.edges(1).boundaryType == 'P') && ...
					~(self.edges(2).boundaryType == 'P')
				
				% store free nodes
				P.free.edge = self.edges(1).nodes(2:end-1);

				% store corresponding replica nodes
				P.replica.edge = flip(self.edges(3).nodes(2:end-1));

			end

			% if East-West periodicity
			if ~(self.edges(1).boundaryType == 'P') && ...
					(self.edges(2).boundaryType == 'P')
				
				% store free nodes
				P.free.edge = self.edges(2).nodes(2:end-1);

				% store corresponding replica nodes
				P.replica.edge = flip(self.edges(4).nodes(2:end-1));

			end

			% if both periodicities
			if (self.edges(1).boundaryType == 'P') && ...
					(self.edges(2).boundaryType == 'P')
				
				% store free edge node
				Pfree_S  = self.edges(1).nodes(2:end-1);
				Pfree_E  = self.edges(2).nodes(2:end-1);
				P.free.edge = [Pfree_S, Pfree_E];

				% store free corner node
				Pfree_SE = self.edges(1).nodes(end);
				P.free.corner = Pfree_SE;

				% store corresponding replica edge nodes
				Preplica_N = flip(self.edges(3).nodes(2:end-1));
				Preplica_W = flip(self.edges(4).nodes(2:end-1));
				P.replica.edge = [Preplica_N, Preplica_W];

				% store corresponding replica corner nodes
				Preplica_SW = self.edges(2).nodes(end);
				Preplica_NW = self.edges(3).nodes(end);
				Preplica_NE = self.edges(4).nodes(end);
				P.replica.corner = [Preplica_SW, Preplica_NW, Preplica_NE];	

			end			
		end
	end
end
