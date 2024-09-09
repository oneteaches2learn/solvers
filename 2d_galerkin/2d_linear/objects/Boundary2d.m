classdef Boundary2d

	properties
		dl
		yLim
		xLim
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
				self.dl = dl;
				[self.xLim,self.yLim] = self.set_limits;
				self.nEdges = self.set_nEdges;
				self.edges = self.setEdgeGeometry;
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
		function [xLim,yLim] = set_limits(self)

			xLim = [min(self.dl(2,:)), max(self.dl(3,:))];
			yLim = [min(self.dl(4,:)), max(self.dl(5,:))];

		end

		function nEdges = set_nEdges(self)

			nEdges.total = size(self.dl,2);
			nEdges.outer = sum(self.dl(7,:) == 0);
			nEdges.yLines = 0;
			nEdges.inclusions = size(self.dl,2) - nEdges.outer;

		end

		function edges = setEdgeGeometry(self)

			edges(self.nEdges.total) = BoundaryEdge2d;
			for j = 1:self.nEdges.total;

				% create new edge
				vert = reshape(self.dl(2:5,j),2,2);
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
			if length(bcTypes) == self.nEdges.total
				for i = 1:self.nEdges.total
					self.edges(i).boundaryType = bcTypes(i);
				end

			% if one bcType is prescribed for all internal edges
			else
				for i = 1:4
					self.edges(i).boundaryType = bcTypes(i);
				end
				for i = 5:self.nEdges.total
					self.edges(i).boundaryType = bcTypes(5);
					bcTypes(i) = bcTypes(5);
				end
			end

		end

		function self = setBCConditions(self,bcConds)

			for i = 1:self.nEdges.total
				self.edges(i).boundaryCondition = bcConds{i};
			end

		end

		function self = setEdgeNodes(self,mesh)

			for i = 1:self.nEdges.total
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
			for i = 1:self.nEdges.total
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


		% UTILITY FUNCTIONS
		function self = add_y_line(self,varargin)

			% store variables
			dl = self.dl;
			if nargin == 1
				yBar = mean(self.yLim);
			else
				for i = 1:length(varargin)
					yBar(i) = varargin{i};
				end
			end
			
			% split dl into outer edges, yLines, and inclusion edges
			dl_outer = dl(:,1:self.nEdges.outer);
			dl_yLines = dl(:,self.nEdges.outer+1:self.nEdges.outer+self.nEdges.yLines);
			dl_inclusions = dl(:,self.nEdges.outer+self.nEdges.yLines+1:end);

			% handle outer edges
			for j = 1:length(yBar)

				% split outer edges per yBar
				for i = [size(dl_outer,2):-1:1]

					% if y_bar is between y_start and y_stop, then split the region
					if (yBar(j) > dl_outer(4,i) && yBar(j) < dl_outer(5,i)) || ...
					(yBar(j) > dl_outer(5,i) && yBar(j) < dl_outer(4,i))
						
						% duplicate the i-th column
						dl_outer = [dl_outer(:,1:i), dl_outer(:,i:end)];
						
						% adjust y_start and y_stop
						dl_outer(5,i) = yBar(j);
						dl_outer(4,i+1) = yBar(j);
					end

				end

				% increment region ids for those outer edges above yBar
				increment_cols = find(sum(dl_outer(4:5,:) > yBar(j),1) > 0);
				dl_outer(6,increment_cols) = dl_outer(6,increment_cols) + 1;

				% increment region ids for those inclusion edges above yBar
				increment_cols = find(sum(dl_inclusions(4:5,:) > yBar(j),1) > 0);
				dl_inclusions(7,increment_cols) = dl_inclusions(7,increment_cols) + 1;

			end

			% update number of outer edges
			self.nEdges.outer = size(dl_outer,2);
			
			% build matrix for yLines
			yChanges = find(dl_outer(4,1:self.nEdges.outer/2) ~= ...
										dl_outer(5,1:self.nEdges.outer/2));
			yChanges = yChanges(2:end);

			dl_yLines =  ...
				[2 * ones(1,length(yChanges));
				self.xLim(1) * ones(1,length(yChanges));
				self.xLim(2) * ones(1,length(yChanges));
				dl_outer(4,yChanges);
				dl_outer(4,yChanges);
				dl_outer(6,yChanges-1);
				dl_outer(6,yChanges);
				zeros(3,length(yChanges))];
			
			% store result
			self.nEdges.outer = size(dl_outer,2);
			self.nEdges.yLines = size(dl_yLines,2);
			self.nEdges.inclusions = size(dl_inclusions,2);
			self.dl = [dl_outer, dl_yLines, dl_inclusions];

			%self.boundary.edges = self.setEdgeGeometry;
			%}

		end


		% PLOTTING FUNCTIONS
		function h = plot(self,NameValueArgs)

			arguments
				self
				NameValueArgs.EdgeLabels string = "off"
				NameValueArgs.FaceLabels string = "off"
			end
			x = NameValueArgs;

			h = pdegplot(self.dl,FaceLabels=x.FaceLabels,EdgeLabels=x.EdgeLabels)

		end

	end
end
