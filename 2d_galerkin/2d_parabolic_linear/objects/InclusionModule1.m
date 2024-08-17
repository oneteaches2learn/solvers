classdef InclusionModule1 < InclusionModule
% INCLUSIONMODULE1() creates inclusions where the ratio of |partial Q| / |Y|
% set by the unit cell. 
%
% Author: Tyler Fara				Date: July 3, 2024
%-----------------------------------------------------------------------------%
% Notes
%	(1) Per [Deufelhard,Hochmuth 2004], the bioheat transfer equation is
%	derived by a homogenization procedure in which inclusion Q is embedded in
%	unit cell Y. This unit cell is then scaled and copied to produce a set of
%	inclusions called Q_eps, a Robin boundary condition is placed on the
%	boundary of Q_eps, and the solution of a sequence of such problems
%	converges to a solution of the bioheat transfer equation, where the
%	coefficient scaling the energy exchange term is 
%
%						incRatio = |partial Q| / |Y|
%
%	the purpose of the InclusionModule objects is to dictate the size of the
%	inclusions comprising Q_eps. Various methods will be tested for generating
%	and scaling the inclusions. So various InclusionModule objects will be
%	created to scale the inclusions in different ways. 
%
%	(2) An InclusionModule object will be an input for the PuncturedDomain2d
%	object. The InclusionModule object will have a set of rules for how the
%	inclusions in the PuncturedDomain will be created. The rules in the
%	InclusionModule object might require certain parameters; so an instance of
%	the InclusionModule object can be created with certain parameters set, then
%	that instance will be stored in the PuncturedDomain2d object. The
%	PuncturedDomain2d object will call the function
%	getGeometryDescriptionMatrix, which will create gd, a geometry description
%	matrix, based on the geometry of the overall domain, the number of
%	inclusions, and the parameters stored in the InclusionModule instance. The
%	output will be the geometry description matrix. 
%
%	(3) The InclusionModule1 object sets incRatio = |partial Q| / |Y| per the
%	following procedure:
%
%		(*) Assume that Omega is a rectangular domain with dimensions w x h.
%		(*) Let Nx and Ny be the number of inclusions in the x and y
%		directions.
%		(*) Assume Y has dimensions (w / Nx) x (h / Ny). 
%		(*) Assume Q is rectangular and the height-width ratio of Q matches
%		that of Y.
%		(*) Set dimensions of Q so that |partial Q| / |Y| = incRatio. Don't
%		worry about keeping incRatio constant as you scale. 
%
%	(4) JULY 3: For now, we assume that Omega is square (and therefore so also
%	are Y and Q). We assume that Omega = (0,1) x (0,1).
%-----------------------------------------------------------------------------%

	properties
		incRatio	
	end

	methods
		function self = InclusionModule1(incRatio)

			self.incRatio = incRatio;
		
		end

		function gd = getGeometryDescriptionMatrix(self,xBounds,yBounds,N_x,N_y)

			% initialize storage
			gd = zeros(10,N_x*N_y+1);

			% put boundary rectangle in column 1 of gd
			gd(:,1) = PuncturedDomain2d.boundsToGeometryDescriptionColumn(xBounds,yBounds);

			% if nonzero number of inclusions
			if (N_x ~= 0) && (N_y ~= 0)

				% get coordinates of puncture in unit cell
				q1 = (1/2 - 1/8 * self.incRatio);
				q2 = (1/2 + 1/8 * self.incRatio);

				% scale coordinates
				q1 = q1 / N_x;
				q2 = q2 / N_x;

				% copy coordinates to create coordinate frame
				coord_base = [q1 q2];
				coord = coord_base;
				for i = 1:N_x-1
					coord_i = coord_base + i/N_x;
					coord = [coord, coord_i];
				end
				coord = [coord 1];
				
				% get x and y coordinates of inclusions
				xVert = coord;
				yVert = coord;

				% put inclusions in remaining columns of gd
				for i = 0:N_y-1, for j = 0:N_x-1

					% current column number
					col = i*N_x + j + 2;

					% bounds for inclusion
					xPunc = [xVert(2*j+1) xVert(2*j+2)];
					yPunc = [yVert(2*i+1) yVert(2*i+2)];

					% add column to geometry description matrix
					gd(:,col) = PuncturedDomain2d.boundsToGeometryDescriptionColumn(xPunc,yPunc);

				end, end
			end
		end

	end

end


