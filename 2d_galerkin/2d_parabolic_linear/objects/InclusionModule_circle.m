classdef InclusionModule_circle

	properties
		incRatio
		epsilon
		Y
	end

	properties (Hidden)
		nCopies = 25
		xRows = [2 3 8]
		yRows = [4 5 9]
		rRow  = [10]
	end

	properties (Hidden,SetAccess = private)
		Y_area
		Q_radius
		Q_circumference
		dl_inclusion
		gm_inclusion
	end

	methods
		% CONSTRUCTOR
		function self = InclusionModule_circle(incRatio,epsilon)

			self.incRatio = incRatio;
			self.epsilon = epsilon

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
			dl = repmat(self.dl_inclusion,1,self.nCopies^2)
			dl(self.xRows,:) = dl(self.xRows,:) + x_translate;
			dl(self.yRows,:) = dl(self.yRows,:) + y_translate;

			% scale dl
			dl([self.xRows,self.yRows,self.rRow],:) = ...
					self.epsilon * dl([self.xRows,self.yRows,self.rRow],:);

		end


		% GETTERS
		function rad = get.Q_radius(self)

			rad = self.Q_circumference / 2 / pi;

		end

		function Y_area = get.Y_area(self)

			Y_area = abs(self.Y(1,1) - self.Y(1,2)) * abs(self.Y(2,1) - self.Y(2,2));

		end

		function Q_circumference = get.Q_circumference(self)

			Q_circumference = self.incRatio * self.Y_area;

		end

		function gm_inclusion = get.gm_inclusion(self)

			xCoord = mean(self.Y(1,:));
			yCoord = mean(self.Y(2,:));

			gm_inclusion = [1,xCoord,yCoord,self.Q_radius,0,0,0,0,0,0]';

		end

		function dl_inclusion = get.dl_inclusion(self)

			dl_inclusion = decsg(self.gm_inclusion);

		end

	end
end
