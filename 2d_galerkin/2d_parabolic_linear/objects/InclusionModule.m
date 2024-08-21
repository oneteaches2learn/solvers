classdef (Abstract) InclusionModule

	properties
		incRatio
		epsilon
		Y
		Q
	end

	properties (Hidden)
		nCopies = 25
		Y_area
		Qeps_centers
	end

	methods
		function self = InclusionModule(Y,incRatio,epsilon)

			self.incRatio = incRatio;
			self.epsilon = epsilon;

			self.Y = [0 1; 0 1];

		end

		% UTILITY FUNCTIONS
		function dl = decomposedMatrix(self)

			% get all permutations of x and y translations
			translate = [0:1:self.nCopies-1];
			combos = table2array(combinations(translate,translate))';

			% copy each column of combos four times
			combos = repmat(combos,4,1);
			combos = reshape(combos,2,[]);

			% isolate x- and y-translations
			x_translate = combos(1,:);
			y_translate = combos(2,:);

			% make unscaled dl
			dl = repmat(self.dl_inclusion,1,self.nCopies^2);
			dl(self.xRows,:) = dl(self.xRows,:) + x_translate;
			dl(self.yRows,:) = dl(self.yRows,:) + y_translate;

			% scale dl
			dl([self.xRows,self.yRows,self.rRow],:) = ...
					self.epsilon * dl([self.xRows,self.yRows,self.rRow],:);

		end


		function Y_area = get.Y_area(self)

			Y_area = abs(self.Y(1,1) - self.Y(1,2)) * abs(self.Y(2,1) - self.Y(2,2));

		end

		function Qeps_centers = get.Qeps_centers(self)

			

		end

	end


	methods (Static)
		function gd = gd_from_dl(dl)


		end
	end

end
