classdef Coefficients

	properties
		porosity
	end

	methods
		function self = Coefficients(por)
			self.porosity = por;
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

			usedVars = Coefficients.getFunctionVariables(f);
			result = any(strcmp(usedVars,checkVar));

		end

		function result = isTimeVarying(f)

			result = Coefficients.hasVariable(f,'t');

		end
	end

end

