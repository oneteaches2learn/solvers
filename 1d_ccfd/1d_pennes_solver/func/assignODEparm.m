function [v_0,mu,c_B,S,vLower,vUpper] = assignODEparm(odeParm)

	v_0 = odeParm{1};
	mu  = odeParm{2};
	c_B = odeParm{3};
	S   = odeParm{4};
	vLower = odeParm{5};
	vUpper = odeParm{6};

end
