classdef Coefficients

	methods
		function self = Coefficients
		end
	end

	methods (Static)
		function vars = getFunctionVariables(fh)

			% use regex to collect variable names from function_handle
			code = func2str(fh);
			pattern = '@\((.*?)\)';
			tokens = regexp(code,pattern,'tokens');

			if ~isempty(tokens)
				vars = strsplit(tokens{1}{1},',');
				vars = strtrim(vars);

			else
				vars = {};

			end
		end

		function result = hasVariable(f,checkVar)

			if isa(f,'double')
				x = sym('x',[1 2]);
				f = symfun(f,x);
			end
			
			if isa(f,'sym') || isa(f,'symfun')
				f = matlabFunction(f,x);
			end
			
			usedVars = Coefficients.getFunctionVariables(f);
			result = any(strcmp(usedVars,checkVar));

		end

		function result = isTimeVarying(f)

			if isa(f,'double')
				x = sym('x',[1 2]);
				f = symfun(f,x);
			end
			
			if isa(f,'sym') || isa(f,'symfun')
				f = matlabFunction(f,x);
			end

			result = Coefficients.hasVariable(f,'t');

		end

		function result = isNonlinear(f)

			if isa(f,'double')
				x = sym('x',[1 2]);
				f = symfun(f,x);
			end

			if isa(f,'sym') || isa(f,'symfun')
				f = matlabFunction(f);
			end
			
			result = Coefficients.hasVariable(f,'u');

		end

		function [f,t,u] = checkVariables(f)
		% Ensures that the function handle f has the correct number of variables.

			if ~isNonlinear(f) && ~isTimeVarying(f)
				f = @(x1,x2,t,u)(f(x1,x2));
				t = 0;
				u = 0;
			elseif ~isNonlinear(f) && isTimeVarying(f)
				f = @(x1,x2,t,u)(f(x1,x2,t));
				u = 0;
			elseif isNonlinear(f) && ~isTimeVarying(f)
				f = @(x1,x2,t,u)(f(x1,x2,u));
				t = 0;
			end

		end

	end

end

