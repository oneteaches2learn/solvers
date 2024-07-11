function [endpoints_L,endpoints_R] = set_subdomain_endpoints(a,b,punctures)

	subdom_num = length(punctures) + 1;

	% Set endpoints of each subdomain
	endpoints_L(1) = a;
	for i = 2:subdom_num
		endpoints_R(i-1) = punctures{i-1}.endpoint_L;
		endpoints_L(i)   = punctures{i-1}.endpoint_R;
	end
	endpoints_R(subdom_num) = b;

end
