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
		effectiveRegion = 'Omega_eps'
		nEdges
		edges
		D_nodes
		N_nodes
		R_nodes
		P_nodes
		freeNodes
	end

	properties (Dependent)
%		nOuterSegments
%		nOuterEdges
%		nInclusionSegments
%		nInclusionEdges
%		nyLines
%		dlIDs_outer
%		dlIDs_inclusions
%		dlIDs_yLines
	end

	properties (Dependent,Hidden)
		nEdges_outer
		nEdges_inclusions
		nEdges_inclusionsRAW
%		dl_outer
%		dl_yLines
%		dl_inclusions
		faceIDs
	end

	properties (Hidden)
		dl
		bcTypes
		bcConds
	end

	methods
		% CONSTUCTOR
		function self = Boundary2d(dl)

			% create edge geometry
			if nargin > 0
				self.dl = dlObj(dl);
				self.edges = self.setEdges;
			end

		end


		% GETTERS 
		function val = get.nEdges(self)

			val = self.nEdges_outer + self.nEdges_inclusions;

		end

		function val = get.nEdges_outer(self)

			val = 4;

		end

		function val = get.nEdges_inclusions(self)

			if strcmp(self.effectiveRegion,'Omega_eps')
				val = self.nEdges_inclusionsRAW;
			else
				val = 0;
			end

		end

		function val = get.nEdges_inclusionsRAW(self)

			val = length(self.dl.edgeSegDict) - self.nEdges_outer;

		end

		function val = get.faceIDs(self)

			faceIDs.Omega_eps = self.dl.faceIDs_Omega_eps;
			faceIDs.Q_eps = self.dl.faceIDs_Q_eps;
			faceIDs.Omega = self.dl.faceIDs_all;
			val = faceIDs;

		end


		% SETTERS
		function edges = setEdges(self)

			edges(self.nEdges) = BoundaryEdge2d;
			for j = 1:self.nEdges;

				% create new edge
				vert = self.dl.edgeCoords(:,:,j);
				edges(j).vertex1 = vert(1,:);
				edges(j).vertex2 = vert(2,:);
				edges(j).segIDs = self.dl.edgeSegDict{j};

			end

		end

		function self = setBCTypes(self,bcTypes)

			% store bcTypes as hidden property
			self = self.storeBCTypes(bcTypes);
			bcTypes = self.bcTypes;

			% set bcTypes on edges
			for i = 1:self.nEdges
				self.edges(i).boundaryType = bcTypes(i);
			end

		end

		function self = storeBCTypes(self,bcTypes)

			% convert cell array to char array if necessary
			if isa(bcTypes,'cell')
				for i = length(bcTypes)
					temp(i) = bcTypes{i};
				end
				bcTypes = temp;
			end
			
			% convert strings to char array if necessary
			if isa(bcTypes,'string')
				bcTypes = char(bcTypes);
			end

			% if one bcType is prescribed for all internal edges, repeat it
			if length(bcTypes) - self.nEdges_outer == 1
				bcTypes = [bcTypes(1:4), bcTypes(end)*ones(1,self.nEdges_inclusionsRAW)];	
			end

			% store bcTypes as hidden property
			self.bcTypes = bcTypes;

		end

		function self = setBCConds(self,bcConds)

			% store bcConds as hidden property
			self = self.storeBCConds(bcConds);
			bcConds = self.bcConds;

			% if one bcCond is prescribed for all internal edges
			if (length(bcConds) - self.nEdges_outer == 1) && ...
										(self.nEdges_inclusions ~= 0)

				% set bcConds on outer edges
				for i = 1:self.nEdges_outer
					self.edges(i).boundaryCondition = bcConds{i};
				end

				% set bcConds on inclusion edges
				for i = 1:self.nEdges_inclusions
					self.edges(i+self.nEdges_outer).boundaryCondition = bcConds(end);
				end

			% else, set bcConds on all appropriate
			else
				for i = 1:self.nEdges
					self.edges(i).boundaryCondition = bcConds{i};
				end
			end
		end
		
		function self = storeBCConds(self,bcConds)

			% store bcTypes as hidden property
			self.bcConds = bcConds;

		end


		% MANIPULATE BOUNDARY NODES
		function self = setEdgeNodes(self,mesh)

			for i = 1:self.nEdges
				edgeID = self.edges(i).segIDs;
				self.edges(i).nodes = mesh.Mesh.findNodes('region','Edge',edgeID);
				self.edges(i).nNodes = length(self.edges(i).nodes);
			end

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

		function self = setFreeNodes(self,mesh)

			boundNodes = [self.D_nodes, self.P_nodes.replica.edge, ...
												self.P_nodes.replica.corner];
			self.freeNodes = setdiff(mesh.effectiveNodes,boundNodes);

		end


		% Y-LINE
		function self = add_yline(self,varargin)

            % update dl
            self.dl = self.dl.add_yline(varargin{:});

			% update edge segment IDs
			self = self.set_edgeSegIDs;

		end

		function self = set_edgeSegIDs(self)

			for i = 1:self.nEdges
				self.edges(i).segIDs = self.dl.edgeSegDict{i};
			end

		end


		% INCLUSIONS ON/OFF
		function self = inclusionsON(self)

			% preserve boundary condition data, if necessary
			if ~isempty(self.bcTypes), bcTypes = self.bcTypes;
			else, bcTypes = []; end
			if ~isempty(self.bcConds), bcConds = self.bcConds;
			else, bcConds = []; end
			if ~isempty(self.edges(end).boundaryCondition),
				for i = 1:self.nEdges
					edgeBCConds{i} = self.edges(i).boundaryCondition;
				end
			else, edgeBCConds = []; end

			% Turn inclusions on
			self.dl = self.dl.inclusionsON;
			self = Boundary2d(self.dl.mat);
			self.effectiveRegion = 'Omega';
			self.edges = self.setEdges;

			% reset boundary conditions
			if ~isempty(bcTypes), self = self.setBCTypes(bcTypes); end
			if ~isempty(bcConds), self = self.setBCConds(bcConds); end
			if ~isempty(edgeBCConds)
				for i = 1:self.nEdges
					self.edges(i).boundaryCondition = edgeBCConds{i};
				end
			end
			
		end

		function self = inclusionsOFF(self)

			% preserve boundary condition data, if necessary
			if ~isempty(self.bcTypes), bcTypes = self.bcTypes;
			else, bcTypes = []; end
			if ~isempty(self.bcConds), bcConds = self.bcConds;
			else, bcConds = []; end
			if ~isempty(self.edges(end).boundaryCondition),
				for i = 1:self.nEdges
					edgeBCConds{i} = self.edges(i).boundaryCondition;
				end
			else, edgeBCConds = []; end

			% Turn inclusions on
			self.dl = self.dl.inclusionsON;
			self = Boundary2d(self.dl.mat);
			self.effectiveRegion = 'Omega_eps';
			self.edges = self.setEdges;

			% reset boundary conditions
			if ~isempty(bcTypes), self = self.setBCTypes(bcTypes); end
			if ~isempty(bcConds), self = self.setBCConds(bcConds); end
			if ~isempty(edgeBCConds)
				for i = 1:self.nEdges
					self.edges(i).boundaryCondition = edgeBCConds{i};
				end
			end

		end


		% PLOTTING FUNCTIONS
		function h = plot(self,NameValueArgs)

			arguments
				self
				NameValueArgs.EdgeLabels string = "off"
				NameValueArgs.FaceLabels string = "off"
			end
			x = NameValueArgs;

			h = pdegplot(self.dl.mat,FaceLabels=x.FaceLabels,EdgeLabels=x.EdgeLabels)

		end

	end
end
