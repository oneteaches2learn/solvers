function y = weightFunction(x,wLower,wUpper)
%WEIGHTFUNCTION(X,WLOWER,WUPPER) creates the weight function defined
%				
%				( 0, x < wLower
%		w(x) =  | linear, wLower < x < wUpper
%				( 1, x > wUpper.
%
% Author: Tyler Fara				Date: April 3, 2024
%-----------------------------------------------------------------------------%
% Inputs
%	x			double or vector of doubles
%	wLower		double
%	wUpper		double
%
% Outputs
%	y			double or vector of doubles
%-----------------------------------------------------------------------------%
% Notes
%	(1) It is intended to make specific, anonymous weight functions by picking
%	a specific wLower and wUpper, then using the syntax
%
%			W = @(x)(weightFunction(x,wLower,wUpper));
%
%	Then, such custom weight functions can be passed as function_handles to
%	other functions that are to be weighted and are designed to take such
%	weight functions as parameters.
%-----------------------------------------------------------------------------%

	y = zeros(size(x));

	for i = 1:length(x)
		if x(i) > wUpper
			y(i) = 1;
		elseif x(i) < wLower
			y(i) = 0;
		else
			y(i) = (x(i) - wLower) / (wUpper - wLower);
		end
	end
	
end


