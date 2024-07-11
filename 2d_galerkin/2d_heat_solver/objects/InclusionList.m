classdef InclusionList
% NOTE: NOT CURRENTLY USED (JULY 1 2024) PROBABLY YOU CAN DELETE THIS

	properties
		Inclusions
		nInclusions
	end

	methods
		function self = InclusionList(gd)

			self.nInclusions = size(gd,2) - 1;
			self.Inclusions = self.generateInclusions(gd);
			
		end

		function inclusions = generateInclusions(self,gd);

			inclusions = [];
			for i = 1:self.nInclusions
				inclusion_i = RectangularInclusion2d(gd(:,i+1));
				inclusions = [inclusions inclusion_i];
			end

		end

	end

end
