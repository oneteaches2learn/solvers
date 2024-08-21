classdef InclusionModule_circle < InclusionModule

	properties
	end

	properties (Hidden)
		xRows = [2 3 8]
		yRows = [4 5 9]
		rRow  = [10]
	end

	properties (Hidden)
		Q_radius
		Q_center
		Q_circumference
		dl_inclusion
		gm_inclusion
	end

	methods
		% CONSTRUCTOR
		function self = InclusionModule_circle(incRatio,epsilon)

			self = self@InclusionModule(incRatio,epsilon);

		end

		% GETTERS
		function rad = get.Q_radius(self)

			rad = self.Q_circumference / 2 / pi;

		end

		function Q_circumference = get.Q_circumference(self)

			Q_circumference = self.incRatio * self.Y_area;

		end

		function Q_center = get.Q_center(self)

			Q_center = mean(self.Y,2);
		
		end

		function gm_inclusion = get.gm_inclusion(self)

			gm_inclusion = [1,self.Q_center(1),self.Q_center(2),self.Q_radius,0,0,0,0,0,0]';

		end

		function dl_inclusion = get.dl_inclusion(self)
		% dl_inclusion uses MATLAB's decsg function to convert the inclusion's
		% geometry description column into a decomposed geometry description
		% matrix. Note that rows 6 and 7 describe which side of the inclusion
		% is "inside" versus "outside" the domain. By default, if a shape's
		% boundary is traversed counterclockwise, MATLAB assumes the region to
		% the left of the boundary is the interior. However, for an inclusion,
		% it is actually the region to the right that is the interior of the
		% domain. So, rows 6 and 7 must be interchanged for an inclusion.

			dl_inclusion = decsg(self.gm_inclusion);
			dl_inclusion([6 7],:) = dl_inclusion([7 6],:);

		end

	end
	
	methods (Static)
		function nor = outward_normal(nod1,nod2)

			seg = [nod1 nod2];
			seg = [diff(seg(1,:));diff(seg(2,:))];
			seg = seg / vecnorm(seg);

			nor = [seg(2); -seg(1)];

		end

	end
end
