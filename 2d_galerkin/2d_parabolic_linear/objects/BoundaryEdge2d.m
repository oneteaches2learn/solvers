classdef BoundaryEdge2d

	properties
		ID
		vertex1
		vertex2
		outwardNormal
		boundaryType
		boundaryCondition
		nodes
		nNodes
	end

	methods
		function self = BoundaryEdge2d(vertex1,vertex2)

			self.vertex1 = vertex1;
			self.vertex2 = vertex2;

		end

	end
end
