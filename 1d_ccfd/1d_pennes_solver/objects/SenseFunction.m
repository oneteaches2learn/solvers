classdef SenseFunction

	properties
		vector	
		value
	end

	methods
		function obj = SenseFunction(domain,input,cutoff)
			if nargin > 2
				input = obj.check_inputs(domain,input,cutoff);
				obj.vector = obj.compute_sense_vector(domain,input,cutoff);
			else
				input = obj.check_inputs(domain,input);
				obj.vector = obj.compute_sense_vector(domain,input);
			end
				obj.value = obj.compute_sense_value(input);
		end

		function input = check_inputs(obj,domain,input,varargin)
			err = "Inputs to SenseFunction must be a Domain object, a vector, and an optional double.";
			if ~isa(domain,'Domain') || ~isvector(input)
				error(err)
			end

			if (length(varargin) == 1) && ~isa(varargin{1},'double')
				error(err)
			end

			err = "Cell number in Domain must match length of input vector.";
			if domain.cellNum ~= length(input)
				error(err)
			end

			if ~iscolumn(input)
				input = input';
			end
		end

		function vec = compute_sense_vector(obj,domain,input,varargin)
			M = domain.cellNum; 
			dx = domain.cellWidths;
			
			if length(varargin) > 0, 
				cutoff = varargin{1}; 
			else 
				cutoff = -Inf;
			end
			
			vec = sparse(M,1); 

			width = 0;
			for i = 1:M
				if input(i) > cutoff
					width = width + dx(i);
					vec(i) = dx(i) * obj.H(input(i) - cutoff);
				end
			end
			
			if width == 0
				...
			else
				vec = vec / width;
			end
			%}
		end

		function val = compute_sense_value(obj,input)
			val = dot(obj.vector,input);
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
