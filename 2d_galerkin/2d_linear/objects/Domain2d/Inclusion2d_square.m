classdef Inclusion2d_square < Inclusion2d
% Inclusion2d_square represents a square Q embedded in a unit cell Y

	properties (Dependent)
		center
		dy1_dn
		dy2_dn
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

		function dy1_dn = get.dy1_dn(self)

			dy1_dn = {0,1,0,-1};

		end

		function dy2_dn = get.dy2_dn(self)

			dy2_dn = {-1,0,1,0};

		end
	end
end
