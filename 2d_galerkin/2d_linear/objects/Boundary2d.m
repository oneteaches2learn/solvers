classdef Boundary2d

	properties
		edges
		nEdges
		boundaryTypes
		boundaryConditions
	end

	methods
		function self = Boundary2d(boundaryTypes,boundaryConditions)

			if nargin == 2
				self.boundaryTypes = boundaryTypes;
				self.boundaryConditions = boundaryConditions;
			end

		end

		function bCond = getBoundaryCondition(self,edgeID)
			bCond = self.boundaryConditions{edgeID};
		end

		function BCnodes = sortBoundaryNodes(self,bNodes,nearestEdge)

			BCnodes.D = [];
			BCnodes.N = [];
			BCnodes.R = [];
			for i = 1:size(bNodes,1)
				if self.boundaryTypes{nearestEdge(i)} == 'D'
					BCnodes.D = [BCnodes.D; bNodes(i,:)];
				elseif self.boundaryTypes{nearestEdge(i)} == 'N'
					BCnodes.N = [BCnodes.N; bNodes(i,:)];
				elseif self.boundaryTypes{nearestEdge(i)} == 'R'
					BCnodes.R = [BCnodes.R; bNodes(i,:)];
				end
			end
			BCnodes.D = unique(BCnodes.D);
		end

		function self = addBoundaryConditions(self,bConds)
			self.boundaryConditions = bConds;
		end

	end
end
