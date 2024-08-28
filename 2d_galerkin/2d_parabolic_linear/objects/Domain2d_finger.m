classdef Domain2d_finger < Domain2d_onOff
% Domain2d_onOff is a domain with inclusions that can be turned on or off
%	
% Author: Tyler Fara					Date: August 21, 2024
%-----------------------------------------------------------------------------%

	methods
		% CONSTRUCTOR
		function self = Domain2d_finger(x,y,inc,eps)
		
			% call Domain2d superclass constructor
			self@Domain2d_onOff(x,y,inc,eps);

			% store geometry matrix
			self.dl = Domain2d_finger.finger_modifier(self.dl);

		end
	end

	methods (Static)
		function dl = finger_modifier(dl)

			% extract domain information
			x1 = dl(2,1);
			x2 = dl(2,2);
			y1 = dl(4,2);
			y2 = dl(5,2);
			rad = (y2 - y1) / 2;

			circ_vec = [1,x2,x2,y1,y2,1,0,x2,y1+rad,rad]';

			dl(:,2) = circ_vec;

		end
	end
end
