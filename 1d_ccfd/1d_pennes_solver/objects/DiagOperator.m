classdef DiagOperator < Operator

	methods 
		function obj = DiagOperator(varargin)
			vectors = varargin;
			vectors = obj.check_inputs(vectors);
			obj.matrix = obj.compute_matrix(vectors);
		end

		function vectors = check_inputs(obj,vectors)
			err = "DiagOperator must take at least one input";
			if length(vectors) == 0
				error(err)
			end

			err = "Inputs to DiagOperator must be vectors and have the same length.";
			M = length(vectors{1});
			for i = 1:length(vectors)
				vec = vectors{i};
				if (~isvector(vec)) || (length(vec) ~= M)
					error(err)
				end
				if ~iscolumn(vec)
					vectors{i} = vec';
				end
			end
		end

		function mat = compute_matrix(obj,vectors)
			product = vectors{1};
			M = length(product);
			
			if length(vectors) > 1
				for i = 2:length(vectors)
					vec = vectors{i};
					product = product .* vec;
				end
			end

			mat = spdiags(product,0,M,M);
		end

	end
end
