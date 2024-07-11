function y = integralAverage(u,dx)
	
	y = dot(u,dx) / (sum(dx));

end
