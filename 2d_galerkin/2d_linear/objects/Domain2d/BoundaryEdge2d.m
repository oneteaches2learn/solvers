classdef BoundaryEdge2d

	properties
		segIDs
		vertex1
		vertex2
		outwardNormal
		boundaryType
		boundaryCondition
		boundaryCondition_ddu
		boundaryCondition_corr
		nodes
		nNodes
	end

	methods
		function self = BoundaryEdge2d(vertex1,vertex2)

			if nargin == 2
				self.vertex1 = vertex1;
				self.vertex2 = vertex2;
			end

		end

	end
end
