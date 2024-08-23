classdef Pore2d_circle < Pore2d
% Pore2d_circle represents a circle
%	Pore2d_circle is centered at self.center and has radius self.radius. The
%	decomposed geometry matrix for Pore2d_circle is comprised of four columns,
%	representing four quarters of a full circle. By convention, the first
%	segment starts at the western-most point of the circle. The segments
%	traverse the circle in counterclockwise orientation, each spanning a 90
%	arc. 
%
%	By assumption, Pore2d_circle's intended use case is as an internal
%	boundary inside some larger domain Omega. Therefore, the region to the
%	"right" (as the circle is traversed counterclockwise) is considered
%	"inside" the domain, and the region to the "left" is considered "outside"
%	the domain. Notice that this means that the region that is enclosed by the
%	circle is "outside" of the domain, which is the opposite of one's usual
%	intuition. 
%
%	For this same reason, the unitNormal vector points toward the center of the
%	circle, since this is therefore an "outward normal" relative to the domain
%	in which Pore2d_circle is inserted. 
%------------------------------------------------------------------------------%
%   Decomposed Geometry Matrix Key:
%
%	dl = [   1    1    1    1;   % <~~shape type, 1 == circle 
%		  -rad    0  rad    0;   % <~~starting x coordinate
%	     	 0  rad    0 -rad;   % <~~ending x coordinate
%	 		 0 -rad    0  rad;   % <~~starting y coordinate
%	      -rad    0  rad    0;   % <~~ending y coordinate
%			 0    0    0    0;   % <~~region label left of segment
%			 1    1    1    1;   % <~~region label right of segment
%			 0    0    0    0;   % <~~x coordinate of circle center
%			 0    0    0    0;   % <~~y coordinate of circle center
%		   rad  rad  rad  rad];  % <~~radius
%
%	where rad = self.radius.
%
%	Note that, by assumption, this dl matrix represents a circle that is
%	centered at the origin.
%------------------------------------------------------------------------------%
%	Unit Normal
%	The unitNormal property returns a 2 x 2 matrix. Each column represents the
%	unit normal to one of the semicircles, with the convention:
%
%		Col1: lower semicircle
%		Col2: upper semicircle
%
%	Note that points on each semicircle can be parameterized by the
%	x-coordinate alone. For a unit circle (i.e. the circle centered at the
%	origin with radius = 1), the parameterization is:
%
%		lower semicircle: (x, -sqrt(1 - x^2))
%		upper semicircle: (x,  sqrt(1 - x^2))
%	
%	These are *also* the same formulae for the unit normal vectors that are
%	oriented away from the circle's center. After all: One way to interpret the
%	above parameterization is that the point (x, -sqrt(1 - x^2), for example,
%	represents a vector pointing from the origin to the point (x,-sqrt(1-x^2)).
%	Because this point lies on the unit circle, this vector has unit length.
%	And this vector is colinear to the outward normal at that point. Therefore,
%	this vector simply *is* the unit vector, oriented away from the circle.
%
%	Because Pore2d_circle is intended to be used as an interior boundary to
%	some domain Omega, the vectors that are returned by unitNormal (in the case
%	that the instance of Pore2d_circle has self.radius = 1) are:
%
%		n_lower = -[x; -sqrt(1 - x^2)];
%		n_upper = -[x;  sqrt(1 - x^2)];
%
%	Lastly, in the intended use-case, Pore2d_circle will be embedded in some
%	unit cell Y, which will be copied many times to form a grid of pores inside
%	of Omega, and this grid will be shrunk down by some factor epsilon.
%	Conceptually, Pore2d_circle's center and radius will be set when it is
%	embedded in Y. Then, each copy will be translated to a new center, and also
%	the radius will be rescaled by epsilon. So, the entire set of unitNormals
%	for the whole grid of pores can be obtained by simply composing the
%	unitNormal for the unit circle (recorded above) with the inverse of these
%	translations and scalings. Schematically, here's how it works:
%
%		n_gridPore = n_unitCircle(scale_Y(scale_eps(translate(x)))
%
%	That is, the gridPore is translated back to the origin. The scaling by
%	epsilon is undone. And the scaling to fit inside Y is undone. That results
%	in the unit circle. So if the x coordinate of some point on gridPore can be
%	transformed in these ways, and then that transformed x-coordinate is fed to
%	n_unitCircle, this will return the correct unit vector.
%	
%	In practice, an instance of Pore2d_circle, called Q, is created when an
%	instance of Inclusion2d_circle is created, which stores Y and Q. So Q has
%	radius Q.radius. And scale_Y = 1 / Q.radius. As a result, the property
%	Q.unitNormal has its argument precomposed with scale_Y. I.e., you will see
%	that the argument of Q.unitNormal is already divided by Q.radius. 
%
%	By contrast, it is the PuncturedDomain2d object that has the information
%	about the grid; i.e. PuncturedDomain2d knows epsilon and the coordinates of
%	the center of each scaled-and-translated pore. Given some pore with center
%	(x_c,y_c), we have translate(x) = x - x_c. And scale_eps(x) = x / epsilon.
%	Therefore, scale_eps(translate(x)) = (x - x_c) / epsilon. 
%
%	Long story short: To use Q.unitNormal, presume that x is a point on
%	gridPore. Compute xBar = (x - x_c) / epsilon. Do unitNormal(xBar). 

	properties
		radius
	end

	properties (Dependent)
		circumference
		area
		dl
		unitNormal
	end

	properties (Hidden)
		dl_scaledRows
		dl_xRows
		dl_yRows
		dl_rRows
		unitNormal_lower
		unitNormal_upper
	end


	methods
		% CONSTRUCTOR
		function self = Pore2d_circle(center,radius)

			% call superclass constructor
			self@Pore2d(center);

			% set additional properties
			self.radius = radius;

		end

		% GETTERS
		function area = get.area(self)

			area = pi * self.radius * self.radius;

		end

		function circ = get.circumference(self)

			circ = 2 * pi * self.radius;
			
		end

		function dl = get.dl(self)

			rad = self.radius;

			dl = [   1    1    1    1;   % <~~shape type, 1 == circle 
				  -rad    0  rad    0;   % <~~starting x coordinate
			     	 0  rad    0 -rad;   % <~~ending x coordinate
			 		 0 -rad    0  rad;   % <~~starting y coordinate
			      -rad    0  rad    0;   % <~~ending y coordinate
					 0    0    0    0;   % <~~region label left of segment
					 1    1    1    1;   % <~~region label right of segment
					 0    0    0    0;   % <~~x coordinate of circle center
					 0    0    0    0;   % <~~y coordinate of circle center
				   rad  rad  rad  rad];  % <~~radius

		end

		function output = get.dl_scaledRows(self)

			output = sort([self.dl_xRows self.dl_yRows, self.dl_rRows]);

		end

		function output = get.dl_xRows(self)

			output = [2 3 8];

		end

		function output = get.dl_yRows(self)

			output = [4 5 9];

		end

		function output = get.dl_rRows(self)

			output = 10;

		end 

		function n = get.unitNormal_upper(self)

			x = sym('x',[1 2]);
			n = symfun(-[x(1) / self.radius; sqrt(1-(x(1) / self.radius)^2)],x);

		end

		function n = get.unitNormal_lower(self)

			x = sym('x',[1 2]);
			n = symfun(-[x(1) / self.radius; -sqrt(1-(x(1) / self.radius)^2)],x);

		end

		function n = get.unitNormal(self)

			x = sym('x',[1 2]);
			n = [self.unitNormal_lower, self.unitNormal_upper];

		end
	end

end
