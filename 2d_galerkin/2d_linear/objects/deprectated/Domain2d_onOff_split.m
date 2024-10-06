classdef Domain2d_onOff_split < Domain2d_onOff
% Domain2d_onOff is a domain with inclusions that can be turned on or off
%	
% Author: Tyler Fara					Date: August 21, 2024
%-----------------------------------------------------------------------------%

	methods
		% CONSTRUCTOR
		function self = Domain2d_onOff_split(dom)
		
			% extract necessary parts of input domain
			x = dom.xLim;
			y = dom.yLim;
			inc = dom.inclusion;
			eps = dom.epsilon;
			dl = dom.dl;

			% call Domain2d superclass constructor
			self@Domain2d_onOff(x,y,inc,eps);

			% store decomposed matrix
			self.dl = dl;
			self.edges = dom.edges;

		end

		% GETTERS
		function elemIDs = elements_Omega(self)
			
			mesh = self.Mesh;
			elemIDs_lower = findElements(mesh,"region",Face=1);
			elemIDs_upper = findElements(mesh,"region",Face=2);
			elemIDs = unique([elemIDs_lower,elemIDs_upper]);

		end

		function nodeIDs = nodes_Omega(self)

			mesh = self.Mesh;
			nodeIDs_lower = findNodes(mesh,"region",Face=1);
			nodeIDs_upper = findNodes(mesh,"region",Face=2);
			nodeIDs = unique([nodeIDs_lower,nodeIDs_upper]);

		end

		% UTILITY FUNCTIONS
		function edges = distributeBoundaryNodes(self)

			edges = self.edges;
			for i = 1:length(self.edges)
				edges(i).nodes = self.Mesh.findNodes('region','Edge',self.edgeID(i));
				edges(i).nNodes = length(edges(i).nodes);
			end

		end

		function mesh = generateMesh(self)

			% create decomposed geometry description matrix for split domain
			dl_split = Domain2d_onOff_split.subdivide_dl(self.dl);
			
			% remove duplicate edge
			dl_trim = [dl_split(:,1:4),dl_split(:,6:end)];

			geo = fegeometry(dl_trim);
			geo = geo.generateMesh(Hmax=self.h,GeometricOrder='linear');
			mesh = geo.Mesh;

		end

		function edgeID = edgeID(self,i)

			if i == 1
				edgeID = 1;
			elseif i == 2
				edgeID = [2 5];
			elseif i == 3
				edgeID = 6;
			elseif i == 4
				edgeID = [7 4];
			else
				edgeID = i + 3;
			end

		end

		function dom_on = get_ON(self)

			dom_on = self;
			dom_on.dl = self.dl(:,1:4);
			dom_on.edges = self.edges(1:4);
			dom_on.nEdges = length(dom_on.edges);

		end

		function dom_off = get_OFF(self)

			dom_off = self;
			dom_off.effectiveNodes = dom_off.nodes_Omega;
			dom_off.effectiveElems = dom_off.elements_Omega;
			dom_off.unusedNodes = dom_off.get_unusedNodes;
			dom_off.unusedElems = dom_off.get_unusedElems;
			dom_off.nNodes = length(dom_off.effectiveNodes);
			dom_off.nElems = length(dom_off.effectiveElems);
			dom_off.freeNodes = intersect(dom_off.freeNodes,dom_off.effectiveNodes);

		end

	end

	methods (Static)
		function dl_split = subdivide_dl(dl)

			% get y-midpoint for outer boundary
			y1 = dl(4,2);
			y2 = dl(5,2);
			midpt = y1 + (y2 - y1) / 2;

			% get indices of y1 and y2 in decomposed geometry matrix
			dl_outer = dl(:,1:4);
			y1_ind = find([NaN(3,4); dl_outer(4:5,:); NaN(5,4)] == y1);
			y2_ind = find([NaN(3,4); dl_outer(4:5,:); NaN(5,4)] == y2);

			% make decomposed geometry matrix for lower half
			dl_lower = dl_outer;
			dl_lower(y2_ind) = midpt;
			dl_lower(7,3) = 2.0;

			% make deomposed geometry matrix for upper half
			dl_upper = dl_outer;
			dl_upper(y1_ind) = midpt;
			dl_upper(6,:) = 2.0;
			dl_upper(7,1) = 1.0;

			% manipulate decomposed geometry matrix for inclusions
			nInc = ((size(dl,2)) - 4) / 4;
			if nInc > 0 
				dl_inc = dl(:,5:end);
				dl_inc(7,end/2+1:end) = 2.0;
				dl_inc(6,:) = dl_inc(6,:) + 1;
			else
				dl_inc = [];
			end

			% assemble decomposed geometry matrix for split domain
			dl_split = [dl_lower, dl_upper, dl_inc];

		end


	end
end
