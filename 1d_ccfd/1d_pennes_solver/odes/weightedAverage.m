function y = weightedAverage(u,dx,w)

	if sum(w(u)) == 0
		y = 0;
	else
		y = dot(w(u),u.*dx) / dot(w(u),dx);
	end

end
