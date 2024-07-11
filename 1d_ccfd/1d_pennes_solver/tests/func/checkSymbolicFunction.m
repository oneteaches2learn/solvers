function f = checkSymbolicFunction(f,varargin)
%CHECKSYMBOLICFUNCTION(F,VARARGIN) ensures that input F is a symbolic function
%	with the correct number and type of variables, specified by VARARGIN.
%
% Author: Tyler Fara				Date: April 3, 2024
%-----------------------------------------------------------------------------%
% Inputs:
%	f			double, sym, or symfun; function to be checked
%	varargin	syms; desired variables in final output
%
% Outputs:
%	f			symfun, output converted (if necessary) to symfun object
%-----------------------------------------------------------------------------%
% Examples:
%	(1) Convert sym in one variable to symfun in two variables
%	syms x; syms t; f = 0*x + 1; class(f), 
%	f = checkSymbolicFunction(f,x,t); class(f)
%
%	(2) Convert double to symfun in two variables
%	syms x; syms t; f = 1; class(f)
%	f = checkSymbolicFunction(f,x,t); class(f)
%
%	(3) Preserve original input
%	syms x; syms t; f = 1; class(f)
%	g = checkSymbolicFunction(f,x,t); class(f), class(g)
%
%	(4) Eliminate unnecessary variable
%	syms x; syms t; f(x,t) = x,
%	f = checkSymbolicFunction(f,x)
%-----------------------------------------------------------------------------%
% Notes:
%	(1) You cannot overwrite a double F with a symfun F. Hence, if input F is a
%	double, then F is stored as TEMP, then F is cleared, then a new symfun F(X)
%	is created.  
%-----------------------------------------------------------------------------%

	% Store desired variables as single symbolic array
	for i = 1:length(varargin)
		x(i) = varargin{i};
	end

	% Check type of input f and convert to symfun with desired variables
	if isa(f,'double')
		const = f;
		clear f;
		f(x) = 0*x(1) + const;
	elseif isa(f,'sym') || isa(f,'symfun')
		f(x) = f;
	end

end
