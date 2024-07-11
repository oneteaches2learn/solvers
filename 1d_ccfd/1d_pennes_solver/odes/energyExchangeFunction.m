function y = energyExchangeFunction(u,v,w_u,w_v,C)

	if nargin < 5
		C = 1;
	end

	for i = 1:length(u)
		y(i) = C * w_v(v) * w_u(u(i));
	end

end
