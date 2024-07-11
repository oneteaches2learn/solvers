classdef PuncturedBoundary2d < Boundary2d

	properties
		boundaryTypesInclusions
		boundaryConditionsInclusions
	end

	methods
		function self = PuncturedBoundary2d(bcTypes,bcConds,bcTypesInc,bcCondsInc)

			self@Boundary2d(bcTypes,bcConds);
			self.boundaryTypesInclusions = bcTypesInc;
			self.boundaryConditionsInclusions = bcCondsInc;

		end
	end

end
