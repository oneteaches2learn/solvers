classdef (Abstract) Inclusion2d
% Inclusion 2d represents a pore Q embedded in unit cell Y
%
%	Author: Tyler Fara				Date: 8/20/24
%------------------------------------------------------------------------------%
%	In a homogenization problem, we have a prototypical pore Q embedded in unit
%	cell Y. This is represented by Inclusion2d. The region Y^* = Y \ Q is
%	copied to form an infinite grid, and then scaled by some epsilon. For
%	domain Omega, we have Omega^epsilon = Omega \ epsilon Y^*, where
%	Omega^epsilon is the region on which a PDE is solved, and Q^epsilon = Omega
%	\ epsilon Y^* is the pore space. 
%
%	To implement this procedure, an instance inc of Inclusion2d stores Y and Q,
%	instances of UnitCell2d and Pore2d respectively. Note that Pore2d is an
%	abstract superclass, with subclasses representing specific pore geometries,
%	e.g. a circle, a square, etc. Similarly, Inclusion2d is an abstract
%	superclass, with subclasses having a specific type of Pore2d. 
%
%	This abstraction is necessary because of the property incRatio.
%	Mathematically, 
%
%		incRatio = |partial Q| / |Y|.
%
%	This ratio is useful in homogenization problems. And so, to instantiate
%	Inclusion2d, we directly pass the desired incRatio, along with xLim and
%	yLim. Instances of Inclusion2d use this information to direct the
%	construction of Y and Q appropriately, with xLim and yLim directing the
%	size Y (and therefore also the quantity |Y|, which is the area of Y). Then
%	incRatio, along with |Y|, are used to compute the |partial Q|, which
%	therefore directs the size of Q. Because different Q geometries will have
%	different relationships between |partial Q| and the actual dimensions of Q
%	itself, the construction of Q is left to the subclasses of Inclusion2d.  
%
%	NOTE: As of 8/20/24, the center of Q is presumed to be the center of Y.
%	There is no reason mathematically why Q should be centered at the center of
%	Y; so this is just a convenience until more generality is needed. 
%------------------------------------------------------------------------------%

	properties
		Y
		Q
		incRatio
	end

	properties (Hidden)
		Q_circumference
	end

	properties (Dependent)
		volumeFraction
	end

	methods
		function self = Inclusion2d(xLim,yLim,incRatio)

			% store inputs
			self.incRatio = incRatio;

			% create Y
			self.Y = UnitCell2d(xLim,yLim);

		end

		% GETTERS
		function Q_circ = get.Q_circumference(self)

			Q_circ = self.incRatio * self.Y.area;

		end

		function frac = get.volumeFraction(self)

			measY = self.Y.area;
			measQ = self.Q.area;
			measYstar = measY - measQ;
			frac = measYstar / measY;

		end

	end

end