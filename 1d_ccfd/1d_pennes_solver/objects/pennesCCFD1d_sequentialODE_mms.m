classdef pennesCCFD1d_sequentialODE_mms
% CONVERGENCE STUDY
%
%	Author: Tyler Fara 				Date May 15, 2024

	% BLACK BOX %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% check function variables
	u     = symfun(u,[x t]);
	k     = symfun(k,vars);
	por   = symfun(por,vars);
	gam   = symfun(gam,vars);
	uStar = symfun(uStar,vars);
	alpha_L = symfun(alpha_L,t);
	alpha_R = symfun(alpha_R,t);

	% manufacture BC
	sigma_L = symfun(u(a,t),t);
	sigma_R = symfun(u(b,t),t);
	sigma_L = matlabFunction(sigma_L);
	sigma_R = matlabFunction(sigma_R);
	leftBC  = sigma_L;
	rightBC = sigma_R;

	% manufacture initial condition
	uInit = symfun(u(x,0),x);

	% manufacture source
	u   = symfun(u,vars);
	u_t = symfun(diff(u,t),vars);
	u_x = symfun(diff(u,x),vars);
	q   = symfun(-k*u_x,vars);
	f   = symfun(por * u_t + diff(q,x) + gam * (u - uStar),vars);

	% generate coefficient function_handles
	u     = matlabFunction(u);
	k     = matlabFunction(k);
	por   = matlabFunction(por);
	gam   = matlabFunction(gam);
	uStar = matlabFunction(uStar);
	uInit = matlabFunction(uInit);
	f     = matlabFunction(f);

	% create body/blood temp function
	c_B   = @(v)(C_B);
	R     = @(v,C_B)(resolvent(v,dt,C_B,rLower,rUpper));
	S     = @(u,dx)(weightedAverage(u,dx,w_I));	
	uStar = @(u,v,v_old,dx)(1.0);

	% Assemble Inputs
	parameters = {k,gam,uStar,f,por};
	boundary   = {leftBC,leftType,rightBC,rightType};
	iterParm   = {atol,maxIter,verb};
	odeParm    = {v_0,mu,c_B,0,0};

	for i = 1:4
		M  = base^i;
		nt = (base^(i))^2;
		dt = T / (nt+1);
		domain = {a,b,M};
		time   = {dt,nt,uInit};
		[U{i},xc{i},tGrid{i},V{i}] = pennesCCFD1d_sequentialODE( ...
						domain,parameters,boundary,time,iterParm,odeParm);
	end

	% POST-PROCESS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% optionally, plot result
	if 1 == 0
		% plot initial time step
		ymin = min(U,[],'all');
		ymax = max(U,[],'all');
		ymin = ymin - 0.1*(ymax - ymin);
		ymax = ymax + 0.1*(ymax - ymin);
		plot(xc,U(1,:));
		ylim([ymin ymax]);
		pause();

		% plot remaining time steps
		for i = 2:length(tGrid)
			plot(xc,U(i,:));
			ylim([ymin ymax]);
			pause(0.3);
		end

		pause();
		close();
	end 

	for i = 1:4
		% set up variables
		h  = max(diff(xc{i}));
		dt = max(diff(tGrid{i}));
		[X,T] = meshgrid(xc{i},tGrid{i});
		uTrue = u(X,T);
		Nt = length(tGrid{i});
		N  = length(xc{i});

		% compute error
		err_p = 0;
		for j = 1:Nt, for k = 1:N
			err_p = err_p + h * dt * abs(U{i}(j,k) - uTrue(j,k))^2;
		end, end
		err{i} = sqrt(err_p);

	end

	ratio{1} = [];
	rate{1} = [];
	for i = 1:3
		ratio{i} = err{i} / err{i+1};
		rate{i}  = log(ratio{i}) / log(base);
	end

	rate
end


