classdef UnitCell2d
% UnitCell2d represents a unit cell in a homogenization problem
%
%	Author: Tyler Fara				Date: 8/20/24
%-----------------------------------------------------------------------------%
%	An instance of UnitCell2d is one of two key ingrediates in an instance of
%	Inclusion2d, the other key ingredient being an instance of Pore2d. Call Y
%	the instance of UnitCell2d, Q the instance of Pore2d, and inc the instance
%	of Inclusion2d. Then inc stores Y and Q. Mathematically, we call Y^* = Y \
%	Q. We make an infinite grid of these Y^*, scale by some epsilon, and, for
%	some domain Omega, we do Omega^epsilon = Omega \ epsilon Y^*. Omega^epsilon
%	is porous domain on which some PDE is solved, and Q^epsilon = Omega \
%	Omega^epsilon is the pore space. 
%------------------------------------------------------------------------------%

	properties
		xLim
		yLim
	end

	properties (Dependent)
		area
		diameter
		circumference
		xWidth
		yWidth
		center
	end

	methods
		% CONSTRUCTOR
		function self = UnitCell2d(xLim,yLim)

			self.xLim = xLim;
			self.yLim = yLim;

		end


		% GETTERS
		function width = get.xWidth(self)

			width = diff(self.xLim);

		end

		function width = get.yWidth(self)

			width = diff(self.yLim);

		end

		function area = get.area(self)

			area = self.xWidth * self.yWidth;

		end

		function diam = get.diameter(self)

			self.diameter = max(self.xWidth,self.yWidth);

		end

		function circ = get.circumference(self)

			self.circumference = 2 * self.xWidth + 2 * self.yWidth;

		end

		function center = get.center(self)

			center = sum([self.xLim; self.yLim],2) / 2; 
			center = center';

		end


		% UTILITY FUNCTIONS
		% DEPRECATED: this method is now the responsibility of PuncturedDomain2d
		function Yeps_centers = Yeps_centers(self,epsilon,xCopies,yCopies)

			xCoord = epsilon * ([0:1:xCopies-1] + 0.5 * self.xWidth);
			yCoord = epsilon * ([0:1:yCopies-1] + 0.5 * self.yWidth);

			Yeps_centers = table2array(combinations(yCoord,xCoord));
			Yeps_centers = Yeps_centers(:,[2 1]);
			
		end
		%----------------------------------------------------------------------%


	end
end
