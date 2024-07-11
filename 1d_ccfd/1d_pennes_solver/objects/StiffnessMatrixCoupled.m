classdef StiffnessMatrixCoupled < Operator

	methods
		function obj = StiffnessMatrixCoupled(stiffness,sense,blood,alpha)
			obj.check_inputs(stiffness,sense,blood,alpha);
			obj.matrix = obj.calculate_stiffness_matrix(stiffness,sense,blood,alpha);
		end

		function check_inputs(obj,stiffness,sense,blood,alpha)
			err = "Inputs must be StiffnessMatrixTime, SenseFunction, ColumnVector, and double";
			if ~isa(stiffness,'StiffnessMatrixTime') || ~isa(sense,'SenseFunction') ...
					|| ~isa(blood,'ColumnVector') || ~isa(alpha,'double')
				error(err)
			end
		end

		function mat = calculate_stiffness_matrix(obj,stiffness,sense,blood,alpha)
			[M,temp] = size(stiffness.matrix);
			mat = sparse(M,M);

			mat(1:M,1:M) = stiffness.matrix;
			mat(M+1,1:M) = -alpha * sense.vector;
			mat(1:M,M+1) = -blood.vector';
			mat(M+1,M+1) = alpha + 1;
		end
	end
end
