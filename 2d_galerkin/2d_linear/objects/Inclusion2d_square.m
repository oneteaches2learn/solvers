classdef Inclusion2d_square < Inclusion2d
% Inclusion2d_square represents a square Q embedded in a unit cell Y

	properties (Dependent)
		center
	end

	properties (Hidden)
		Q_sideLength
	end

	methods
		function self = Inclusion2d_square(xLim,yLim,incRatio)

			% call superclass constructor
			self@Inclusion2d(xLim,yLim,incRatio);

			% create Q
			self.Q = Pore2d_square(self.center,self.Q_sideLength);

		end

		function center = get.center(self)

			center = self.Y.center;

		end

		function Q_sideLength = get.Q_sideLength(self)

			Q_sideLength = self.Q_circumference / 4;

		end

	end
end
