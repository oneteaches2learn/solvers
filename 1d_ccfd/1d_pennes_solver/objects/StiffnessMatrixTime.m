classdef StiffnessMatrixTime < Operator

	methods 
		function obj = StiffnessMatrixTime(domain,time,varargin)
			operators = varargin;
			obj.check_inputs(domain,time,operators);
			obj.matrix = obj.assemble_stiffness_matrix(domain,time,operators);
		end
		
		function check_inputs(obj,domain,time,operators)
			err = sprintf(['Inputs to StiffnessMatrixTime(domain,time,varargin) must be ' ...
				   'a Domain object, a TimeStepping object, and optional Operator objects']);
			if ~isa(domain,'Domain') || ~isa(time,'TimeStepping')
				error(err)
			end
			
			if length(operators) > 0
				for i = 1:length(operators)
					if ~isa(operators{i},'Operator')
						error(err)
					end
				end
			end

			err = "Size of operators does not match";
			if length(operators) > 0
				M = domain.cellNum;
				for i = 1:length(operators)
					[height, width] = size(operators{i}.matrix);
					if height ~= M
						error(err)
					end
				end
			end

		end

		function mat = assemble_stiffness_matrix(obj,domain,time,operators)
			M = domain.cellNum;
			tau = time.stepSize;
			mat = spdiags(domain.cellWidths,0,M,M);

			if length(operators) > 0
				for i = 1:length(operators)
					mat = mat + tau * operators{i}.matrix;
				end
			end
		end

	end
end
