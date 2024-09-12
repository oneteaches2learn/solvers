classdef Inclusion2d_circle < Inclusion2d
% Inclusion2d_circle represents a circle Q embedded in unit cell Y

	properties (Dependent)
		center
		unitNormal
		dy1_dn
		dy2_dn
	end

	properties (Hidden)
		Q_radius
		Q_dl
	end

	properties (Hidden,Dependent)
		unitNormal_lower
		unitNormal_upper
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

		function n = get.unitNormal(self,dom)

			n_lower = self.unitNormal_lower;
			n_upper = self.unitNormal_upper;

			n = {n_lower,n_lower,n_upper,n_upper};

		end

		function n = get.unitNormal_lower(self)

			% setup symbolic functions for circle edges
			x = sym('x',[1 2],'real');

			% if unit vectors to circle are needed, store in advance
			n = self.Q.unitNormal_lower;

			% set x-translation
			x_translate = self.center(1);
			x_transformed = (x(1) - x_translate);

			% get normal vector
			n = symfun(n(x_transformed,x(2)),x);

		end

		function n = get.unitNormal_upper(self,dom)

			% setup symbolic functions for circle edges
			x = sym('x',[1 2],'real');

			% if unit vectors to circle are needed, store in advance
			n = self.Q.unitNormal_upper;

			% set x-translation
			x_translate = self.center(1);
			x_transformed = (x(1) - x_translate);

			% get normal vector
			n = symfun(n(x_transformed,x(2)),x);

		end
		
		function dy1_dn = get.dy1_dn(self)

			x = sym('x',[1 2],'real');
			y_1 = symfun([1;0],x);

			n = self.unitNormal;
			for i = 1:4
				dy1_dn{i} = symfun(sum(n{i}.*y_1),x);
			end

		end

		function dy2_dn = get.dy2_dn(self)

			x = sym('x',[1 2],'real');
			y_2 = symfun([0;1],x);

			n = self.unitNormal;
			for i = 1:4
				dy2_dn{i} = symfun(sum(n{i}.*y_2),x);
			end

		end

	end
end
