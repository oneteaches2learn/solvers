function [a,b] = unpack_mms_domain(domain)
	if length(domain) == 2
		a = 0;
		b = domain{1};
	elseif length(domain) == 3
		a = domain{1}; 
		b = domain{2};
	end
end
