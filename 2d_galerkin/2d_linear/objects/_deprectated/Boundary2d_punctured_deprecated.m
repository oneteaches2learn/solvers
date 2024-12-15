classdef Boundary2d_punctured < Boundary2d

	properties
		boundaryTypesInclusions
		boundaryConditionsInclusions
	end

	methods
		function self = Boundary2d_punctured(bcTypes,bcConds,bcTypesInc,bcCondsInc)

			self@Boundary2d(bcTypes,bcConds);
			self.boundaryTypesInclusions = bcTypesInc;
			self.boundaryConditionsInclusions = bcCondsInc;

		end
	end

end
