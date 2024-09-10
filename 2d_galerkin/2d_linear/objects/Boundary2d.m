classdef Boundary2d
% NOTE: The way Boundary2d and Mesh2d interact reflects the state of two
% variables: meshInclusions and effectiveRegion, which are properties of
% Domain2d_punctured. If meshInclusions is 'on', then the interior of inclusions
% are meshed. If effectiveRegion is 'Omega', then the interior of the
% inclusions is included in the effective region. 
%
% In more detail, there are four cases to consider:
%
%  1. The domain has no inclusions. If this is the desired case, then construct
%  a Domain2d object. When using the Domain2d object, the state of the
%  'meshInclusions' and 'effectiveRegion' properties is irrelevant (and indeed
%  these variables are never stored even if, for example, they are included in
%  the mmsParams variable during an mmsTest). In this case, the Boundary2d
%  object will have no inclusions and the effective region will be the entire
%  domain.
%  
%  2. The domain has inclusions and they are not meshed. If this is the desired
%  case, then construct a Domain2d_punctured object. The default state of the
%  'meshInclusions' property is 'off', so this object can be constructed with no
%  additional name-value pair parameters. If one does pass name-value pair
%  parameters to the Domain2d_punctured object, then passing
%  meshInclusions='off' will have the same result. Because the interior of the
%  inclusions is not meshed, the effective region will automatically be
%  Omega_eps, i.e. Omega minus the interior of the inclusions. Thus,
%  effectiveRegion='Omega_eps' is the default state of the 'effectiveRegion'
%  property. However, even if effectiveRegion='Omega' is passed as a name-value
%  pair parameter, this will have no effect. 
%
%  3. The domain has inclusions and they are meshed, but not included in the
%  effective region. This state can be acheived by passing meshInclusions='on'
%  and effectiveRegion='Omega_eps'. Importantly, the boundaries of the inclusion
%  are considered as part of the 'boundary' of the domain, and boundary
%  conditions also need to be set / calculated on these boundary edges. 
%
%  4. The domain has inclusions and they are meshed and included in the
%  effective region. This state can be acheived by passing meshInclusions='on'
%  and effectiveRegion='Omega'. In this case, the boundaries of the inclusion
%  are not considered as part of the 'boundary' of the domain, and boundary
%  conditions do not need to be set / calculated on these boundary edges. Note
%  that in this case, if boundary conditions *are* calculated and stored on the
%  corresponding boundary edges, they will nevertheless be ignored as the nEdges
%  property will be set to 4, and so any loops that try to apply boundary
%  conditions will terminate before reaching the edges that correspond to the
%  inclusions. 
%
% The reason for all of this trouble is that I want to compare the results of a
% homogenized PDE on all of Omega versus that of the sequence of individual PDEs
% just on Omega_eps. Thus, in each such trial, I want to ensure that the mesh in
% Omega_eps is the same whether the inclusions are included in the mesh or not.
% Thus, a mesh is generated over all of Omega, including within the inclusions.
% Then, the effective region is set to Omega_eps, and the inclusions are either
% included in the effective region or not.

	properties
		dl
		dl_IDs
		xLim
		yLim
		nEdges
		edges
		D_nodes
		N_nodes
		R_nodes
		P_nodes
		freeNodes
	end

	properties (Dependent,Hidden)
		dl_outer
		dl_yLines
		dl_inclusions
		faceIDs
	end

	methods
		% CONSTUCTOR
		function self = Boundary2d(dl,bcTypes,bcConds)

			% create edge geometry
			if nargin > 0
				self.dl = self.set_dl(dl);
				[self.xLim,self.yLim] = self.set_limits;
				self.dl_IDs = self.set_dl_IDs;
				self.nEdges = self.set_nEdges;
				self.edges = self.setEdgeGeometry;
			end

			% set boundary condition types
			if nargin > 1 && ~strcmp(bcTypes,'skip')
				self = self.setBCTypes(bcTypes);
			end

			% set boundary conditions
			if nargin > 2 && ~strcmp(bcConds,'skip')
				self = self.setBCConditions(bcConds);
			end

		end

		% SETTERS
		function dl = set_dl(self,dl)

			if size(dl,1) ~= 10
				pad = 10 - size(dl,1);
				dl = [dl; zeros(pad,size(dl,2))];
			end

		end

		function [xLim,yLim] = set_limits(self)

			xLim = [min(self.dl(2,:)), max(self.dl(3,:))];
			yLim = [min(self.dl(4,:)), max(self.dl(5,:))];

		end

		function nEdges = set_nEdges(self,meshInclusions,effectiveRegion)

			if nargin == 1
				nEdges = size(self.dl_outer,2) + size(self.dl_inclusions,2);
			elseif meshInclusions == 'off'
				nEdges = size(self.dl_outer,2) + size(self.dl_inclusions,2);
			elseif meshInclusions == 'on' && effectiveRegion == 'Omega_eps'
				nEdges = size(self.dl_outer,2) + size(self.dl_inclusions,2);
			else
				nEdges = size(self.dl_outer,2);
			end

		end

		function dl_IDs = set_dl_IDs(self)
		% dl_IDs identifies the columns of dl matrix
		%   dl_IDs should be called when the domain object (and therefore also
		%   its boundary) is first instantiated. After creating the domain and
		%   before meshing the domain, additional dividing lines can be added to
		%   the domain which will not be included in the boundary. This is done
		%   using the dom.add_y_line or boundary.add_y_line methods. This method
		%   will update the dl_IDs property.

			dl_IDs.outer = find(self.dl(7,:) == 0);
			dl_IDs.yLines = [];
			dl_IDs.inclusions = setdiff(1:size(self.dl,2),dl_IDs.outer);

		end

		function edges = setEdgeGeometry(self)

			edges(self.nEdges) = BoundaryEdge2d;
			for j = 1:self.nEdges;

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
				edgeID = self.edges(i).ID
				self.edges(i).nodes = mesh.Mesh.findNodes('region','Edge',edgeID);
				self.edges(i).nNodes = length(self.edges(i).nodes);
			end

		end

		function self = setFreeNodes(self,mesh)

			boundNodes = [self.D_nodes, self.P_nodes.replica.edge, ...
												self.P_nodes.replica.corner];
			self.freeNodes = setdiff(mesh.effectiveNodes,boundNodes);

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


		% GETTERS
		function val = get.dl_outer(self)

			val = self.dl(:,self.dl_IDs.outer);

		end

		function val = get.dl_yLines(self)

			val = self.dl(:,self.dl_IDs.yLines);

		end

		function val = get.dl_inclusions(self)

			val = self.dl(:,self.dl_IDs.inclusions);

		end

		function val = get.faceIDs(self)

			faceIDs.Omega_eps = unique(self.dl_outer(6,:));
			faceIDs.Q_eps = unique(self.dl_inclusions(6,:));
			if faceIDs.Q_eps == 0, faceIDs.Q_eps = []; end
			faceIDs.Omega = [faceIDs.Omega_eps, faceIDs.Q_eps];
			val = faceIDs;

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
			dl_outer = self.dl_outer;
			dl_yLines = self.dl_yLines;
			dl_inclusions = self.dl_inclusions;

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

			% build matrix for yLines
			yChanges = find(dl_outer(4,1:end/2) ~= ...
										dl_outer(5,1:end/2));
			yChanges = yChanges(2:end);

			dl_yLines =  ...
				[2 * ones(1,length(yChanges));				% line type
				self.xLim(1) * ones(1,length(yChanges));	% x_start
				self.xLim(2) * ones(1,length(yChanges));	% x_stop
				dl_outer(4,yChanges);						% y_start
				dl_outer(4,yChanges);						% y_stop
				dl_outer(6,yChanges);						% left region ID
				dl_outer(6,yChanges-1);						% right region ID
				zeros(3,length(yChanges))];			 	    % padding

			% update outer edge IDs
			n_outerSegments = size(dl_outer,2)
			outerIDs = 1:n_outerSegments;
			outerIDs = reshape(outerIDs,[],2)';

			southIDs = outerIDs(1,1);
			eastIDs  = outerIDs(1,2:end);
			northIDs = outerIDs(2,1);
			westIDs  = outerIDs(2,2:end);

			self.edges(1).ID = southIDs;
			self.edges(2).ID = eastIDs;
			self.edges(3).ID = northIDs;
			self.edges(4).ID = westIDs;

			% update inclusion IDs
			if length(self.edges) > 4
				for i = 5:length(self.edges)
					j = i - 4;
					self.edges(i).ID = j + n_outerSegments;
				end
			end
			
			% store result
			self.dl = [dl_outer, dl_yLines, dl_inclusions];

			% update dl_IDs
			self.dl_IDs.outer = 1:size(dl_outer,2);
			self.dl_IDs.inclusions = [1:size(dl_inclusions,2)] + (size(self.dl,2) - size(dl_inclusions,2));
			self.dl_IDs.yLines = setdiff([1:size(self.dl,2)],self.dl_IDs.outer);
			self.dl_IDs.yLines = setdiff(self.dl_IDs.yLines,self.dl_IDs.inclusions);

		end

		function self = meshInclusions(self,NameValueArgs)

			arguments
				self
				NameValueArgs.meshInclusions = 'off'
			end

			% store variables
			dl_outer = self.dl_outer;
			dl_yLines = self.dl_yLines;
			dl_inclusions = self.dl_inclusions;

			% get number of faces
			nFaces = max(dl_outer(6:7,:),[],'all');

			% get inclusion numbers to be added
			incFaces = [1:size(dl_inclusions,2)/4] + nFaces;
			incFaces = repmat(incFaces,4,1);
			incFaces = reshape(incFaces,1,[]);

			% optionally, the mesh inclusions can be removed from meshed region
			if strcmp(NameValueArgs.meshInclusions,'off')
				incFaces = incFaces * 0;	
			end

			% add (or remove) inclusions to meshed region
			dl_inclusions(6,:) = incFaces;

			% store result
			self.dl = [dl_outer, dl_yLines, dl_inclusions];

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
