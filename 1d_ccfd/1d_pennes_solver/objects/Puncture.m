classdef Puncture

	properties
		endpoint_L
		endpoint_R
		BC_L
		BC_R
		BCtype_L
		BCtype_R
	end

	methods
		function obj = Puncture(endpoint_L,endpoint_R,BC_L,BC_R,BCtype_L,BCtype_R)
			obj.check_inputs(endpoint_L,endpoint_R,BC_L,BC_R,BCtype_L,BCtype_R)
			
			obj.endpoint_L = endpoint_L;
			obj.endpoint_R = endpoint_R;
			obj.BC_L = BC_L;
			obj.BC_R = BC_R;
			obj.BCtype_L = BCtype_L;
			obj.BCtype_R = BCtype_R;

		end

		function check_inputs(obj,endpoint_L,endpoint_R,BC_L,BC_R,BCtype_L,BCtype_R)
			err1 = "endpoint_L, endpoint_R should be doubles. ";
			err2 = "BC_L, BC_R should be function_handles. ";
			err3 = "BCtype_L,BCtype_R should be chars 'D', 'N', or 'R'.";
			err = strcat(err1,err2,err3);

			if ~isa(endpoint_L,'double') || ~isa(endpoint_R,'double') 
				error(err)
			end
			if (BCtype_L == 'D') || (BCtype_L == 'N')
				if ~isa(BC_L,'function_handle') 
					error(err)
				end
			end
			if (BCtype_R == 'D') || (BCtype_R == 'N')
				if ~isa(BC_R,'function_handle') 
					error(err)
				end
			end
			if ~isa(BCtype_L,'char') || ~isa(BCtype_R,'char') 
				error(err)
			end

			if (sum(BCtype_L == 'DNR') ~= 1) || (sum(BCtype_R == 'DNR') ~= 1)
				error(err)
			end
		end

	end
end
