classdef (Abstract) Pore2d
% Pore2d is an abstract class representing some pore in a porous domain
%
%	Author: Tyler Fara			Date: 8/20/2024
%------------------------------------------------------------------------------%
%	Subclasses of Pore2d represent various geometries of pores, e.g. circles,
%	squares, etc. By Q we mean the region enclosed by the pore. By Y we mean
%	some rectangular unit cell in which Q is embedded. We say that Y^* = Y \ Q.
%	Assume we have some domain Omega. Theoretically, we can make an infinite
%	grid of the Y^*, then scale this set by epsilon. We then intersect Omega
%	and epsilon Y^*, and call the result Omega^epsilon. The pore space, i.e.
%	the part of Omega that has been excluded by this procedure, is called
%	Q^epsilon. In a homogenization problem, we solve a PDE on Omega^epsilon,
%	subject to some boundary conditions on partial Q^epsilon (and partial
%	Omega). 
%
%	To implement this procedure, we introduce the Pore2d class and the
%	UnitCell2d class. By Q we designate an instane of the Pore2d class, and by
%	Y we designate an instance of UnitCell2d. We also introduce the Inclusion2d
%	class, which stores an instance of Y and Q. By inc we mean an instance of
%	Inclusion2d. An instance dom of the PuncturedDomain2d class stores inc,
%	which represents a prototypical pore. dom also stores epsilon and
%	information about Omega. In this way, dom can direct the construction of
%	Q^epsilon, and therefore can direct the construction of Omega^epsilon. 
%
%	Subclasses of Pore2d also have information about the unit normal to Pore2d.
%	By composing the formula for the unit normal with information about the
%	spacing and scaling of pores in Q^epsilon, dom can recover a set of unit
%	normal vectors for the pores in Q^epsilon. 
%------------------------------------------------------------------------------%

	properties
		center
	end

	properties (Abstract)
		circumference
		area
		dl
		unitNormal
	end

	methods
		function self = Pore2d(center)

			self.center = center;

		end
	end

end

