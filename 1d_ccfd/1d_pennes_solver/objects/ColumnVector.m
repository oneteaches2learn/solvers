classdef ColumnVector 

	properties
		vector	
	end

	methods 
		function obj = ColumnVector(varargin)
			vectors = varargin;
			vectors = obj.check_inputs(vectors);
			obj.vector = obj.compute_product(vectors);
		end

		function vectors = check_inputs(obj,vectors)
			err = "ColumnVector must take at least one input";
			if length(vectors) == 0
				error(err)
			end

			err = "Inputs to ColumnVector must be vectors and have the same length.";
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

		function vec = compute_product(obj,vectors)
			product = vectors{1};
			M = length(product);
			
			if length(vectors) > 1
				for i = 2:length(vectors)
					vec = vectors{i};
					product = product .* vec;
				end
			end

			vec = product;
		end

	end
end
