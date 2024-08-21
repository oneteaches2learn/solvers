classdef Inclusion2d_circle < Inclusion2d
% Inclusion2d_circle represents a circle Q embedded in unit cell Y

	properties (Dependent)
		center
	end

	properties (Hidden)
		Q_radius
		Q_dl
	end

	methods
		function self = Inclusion2d_circle(xLim,yLim,incRatio)

			% call superclass constructor
			self@Inclusion2d(xLim,yLim,incRatio);

			% create Q
			self.Q = Pore2d_circle(self.center,self.Q_radius);

		end

		function center = get.center(self)

			center = self.Y.center;

		end

		function Q_radius = get.Q_radius(self)

			Q_radius = self.Q_circumference / 2 / pi;

		end

	end
end
