classdef DiscreteLaplacian < Operator

	methods
		function obj = DiscreteLaplacian(transmissibility)
			obj.check_inputs(transmissibility);
			obj.matrix = obj.build_base_matrix(transmissibility);
		end

		function check_inputs(obj,transmissibility)
			err = "Input for DiscreteLaplacian must be a Transmissibility object.";
			if ~isa(transmissibility,'Transmissibility') 
				error(err)
			end
		end

		function A = build_base_matrix(obj,transmissibility)
		%BUILDBASEMATRIX builds the tridiagonal matrix using transmissibility. 

			tx = transmissibility.cellEdges;
	    	n = length(tx)-1;
	    	A = sparse(n,n);
    
	    	for i = 2:n
	        	gl = i-1;
	        	gr = i;
    
	        	A(gl,gl) = A(gl,gl) + tx(i);
	        	A(gl,gr) = A(gl,gr) - tx(i); 
	        	A(gr,gl) = A(gr,gl) - tx(i);
	        	A(gr,gr) = A(gr,gr) + tx(i);
	    	end
		end

		function print(obj,n)
			if nargin == 1
				full(obj.matrix)
			else	
				mat = obj.matrix;
				full(mat(1:n,1:n))
			end
		end

	end
end
