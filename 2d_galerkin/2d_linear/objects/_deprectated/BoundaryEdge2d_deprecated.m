classdef BoundaryEdge2d_deprecated

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
		function self = BoundaryEdge2d_deprecated(vertex1,vertex2,n)

			if nargin == 2
				self.vertex1 = vertex1;
				self.vertex2 = vertex2;
			end

			if nargin == 3
				self.outwardNormal = n;
			end

		end

	end
end
