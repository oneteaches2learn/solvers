classdef BoundaryCondition

	properties
		condition
		side
		valueLHS
		valueRHS
	end

	methods
		function obj = BoundaryCondition(condition,side)
			side = obj.check_inputs(condition,side);
			obj.side = side;
			obj.condition = condition;
		end

		function side = check_inputs(obj,condition,side)
			err = ['DirichletBoundary inputs must be a double representing the ' ...
					'desired condition, a string "L" or "R", and a ' ...
					'Transmissibility object'];
			
			if ~isa(condition,'double')
				error(err)
			end

			if ~isa(side,'char') && ~isa(side,'string')
				error(err)
			end

			if isa(side,'char')
				side = string(side)
			end

			if ~strcmp(side,"L") && ~strcmp(side,"R")
				error(err)
			end

		end

	end
end
