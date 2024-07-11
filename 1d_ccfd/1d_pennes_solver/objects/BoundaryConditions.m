classdef BoundaryConditions < Operator

	properties
		left
		right
		source
	end

	methods
		function obj = BoundaryConditions(left,right,domain)
			obj.check_inputs(left,right,domain);
			obj.left = left;
			obj.right = right;
			obj.matrix = obj.compute_boundary_matrix(domain);
			obj.source = obj.compute_boundary_source(domain);
		end

		function check_inputs(obj,left,right,domain);
			err = "Inputs to BoundaryConditions must be BoundaryCondition objects and a Domain object.";
			if ~isa(left,'BoundaryCondition') || ~isa(right,'BoundaryCondition') || ~isa(domain,'Domain')
				error(err)
			end

			if ~strcmp(left.side,"L")
				error "Left BoundaryCondition side property is not 'L'";
			end
			if ~strcmp(right.side,"R")
				error "Right BoundaryCondition side property is not 'R'";
			end
		end

		function mat = compute_boundary_matrix(obj,domain)
			n = domain.cellNum;
			mat = sparse(n,n);
			mat(1,1) = obj.left.valueLHS;
			mat(n,n) = obj.right.valueLHS;
		end

		function vec = compute_boundary_source(obj,domain)
			n = domain.cellNum;
			vec = sparse(n,1);
			vec(1) = obj.left.valueRHS;
			vec(n) = obj.right.valueRHS;
		end


		function print(obj)
			matrix = full(obj.matrix)
			source = full(obj.source)
		end
			
	end
end

