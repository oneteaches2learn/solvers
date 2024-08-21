classdef Pore2d_square < Pore2d
% Pore2d_square represents a square
%
%	Author: Tyler Fara			Date: 8/20/2024
%------------------------------------------------------------------------------%
%	Pore2d_square is centered at self.center and has side length
%	self.sideLength. The decomposed geometry matrix for Pore2d_circle is
%	comprised of four columns, representing four sides of the square. By
%	convention, the first segment represents the bottom edge of the square. The
%	segments traverse the circle in counterclockwise orientation.
%
%	By assumption, Pore2d_square's intended use case is as an internal
%	boundary inside some larger domain Omega. Therefore, the region to the
%	"right" (as the square is traversed counterclockwise) is considered
%	"inside" the domain, and the region to the "left" is considered "outside"
%	the domain. Notice that this means that the region that is enclosed by the
%	square is "outside" of the domain, which is the opposite of one's usual
%	intuition. 
%
%	For this same reason, the unitNormal vector points toward the center of the
%	square, since this is therefore an "outward normal" relative to the domain
%	in which Pore2d_square is inserted. 
%------------------------------------------------------------------------------%
%   Decomposed Geometry Matrix Key:
%
%	xLim = self.xLim - self.center(1);
%	yLim = self.yLim - self.center(2); 
%
%	dl = [     2       2       2       2;   % <~~shape type, 2 == line
%	      xLim(1) xLim(2) xLim(2) xLim(1);  % <~~starting x coordinate
%	      xLim(2) xLim(2) xLim(1) xLim(1);  % <~~ending x coordinate
%		  yLim(1) yLim(1) yLim(2) yLim(2);  % <~~starting y coordinate
%		  yLim(1) yLim(2) yLim(2) yLim(1);  % <~~ending y coordinate
%			   0       0       0       0;   % <~~region label left of segment
%			   1       1       1       1;   % <~~region label right of segment
%			   0       0       0       0;   % <~~pad row
%			   0       0       0       0;   % <~~pad row
%		       0       0       0       0];  % <~~pad row
%
%	Note that, by assumption, this dl matrix represents a square centered at
%	the origin. 
%------------------------------------------------------------------------------%
%	unitNormal
%	The unitNormal property returns a 2 x 4 matrix. Each column represents the
%	unit normal vector to one of the sides of Pore2d_square, oriented to face
%	toward the center of the square, and following the convention:
%
%		Col1: bottom edge
%		Col2: right edge
%		Col3: top edge
%		Col4: left edge
%-----------------------------------------------------------------------------%

	properties
		sideLength
	end

	properties (Dependent)
		area
		circumference
		xLim
		yLim
		dl
		unitNormal
	end

	properties (Hidden)
		dl_scaledRows
		dl_xRows
		dl_yRows
		dl_rRows
	end


	methods
		% CONSTRUCTOR
		function self = Pore2d_square(center,sideLength)

			% call superclass constructor
			self@Pore2d(center);

			% set additional propertes
			self.sideLength = sideLength;

		end
		
		% GETTERS
		function area = get.area(self)

			self.area = self.sideLength * self.sideLength;

		end

		function circ = get.circumference(self)

			self.circumference = 4 * self.sideLength;

		end

		function xLim = get.xLim(self)

			xLim = self.center(1) + 0.5 * self.sideLength * [-1,1];

		end

		function yLim = get.yLim(self)

			yLim = self.center(2) + 0.5 * self.sideLength * [-1,1];

		end

		function dl = get.dl(self)

			xLim = self.xLim - self.center(1);
			yLim = self.yLim - self.center(2); 

			dl = [     2       2       2       2;   % <~~shape type, 2 == line
			      xLim(1) xLim(2) xLim(2) xLim(1);  % <~~starting x coordinate
			      xLim(2) xLim(2) xLim(1) xLim(1);  % <~~ending x coordinate
				  yLim(1) yLim(1) yLim(2) yLim(2);  % <~~starting y coordinate
				  yLim(1) yLim(2) yLim(2) yLim(1);  % <~~ending y coordinate
					   0       0       0       0;   % <~~region label left of segment
					   1       1       1       1;   % <~~region label right of segment
					   0       0       0       0;   % <~~pad row
					   0       0       0       0;   % <~~pad row
				       0       0       0       0];  % <~~pad row
		end

		function output = get.dl_scaledRows(self)

			output = sort([self.dl_xRows self.dl_yRows]);

		end

		function output = get.dl_xRows(self)

			output = [2 3];

		end

		function output = get.dl_yRows(self)

			output = [4 5];

		end

		function n = get.unitNormal(self)

			n_lower = [0; 1];
			n_right = [-1; 0];
			n_upper = [0; -1];
			n_left  = [1; 0];

			n = [n_lower, n_right, n_upper, n_left];

		end

	end

end
