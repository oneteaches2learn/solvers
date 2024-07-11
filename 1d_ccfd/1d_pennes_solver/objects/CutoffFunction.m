classdef CutoffFunction

	properties
		vector	
	end

	methods
		function obj = CutoffFunction(input,cutoff)
			obj.check_inputs(input,cutoff);
			obj.vector = obj.compute_cutoff_vector(input,cutoff);
		end

		function input = check_inputs(obj,input,cutoff)
			err = "Inputs to CutoffFunction must be a vector and a double.";
			if ~isvector(input) || ~isa(cutoff,'double')
				error(err)
			end

			if ~iscolumn(input)
				input = input';
			end
		end

		function vec = compute_cutoff_vector(obj,input,cutoff)
			M = length(input);
			vec = sparse(M,1); 

			for i = 1:M
				vec(i) = obj.H(input(i) - cutoff);
			end
		end

		function y = H(obj,x);
			if x > 0
				y = 1;
			else 
				y = 0;
			end
		end
	end
end
