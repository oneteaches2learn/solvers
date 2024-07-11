function v = bodyTempODE(mu,dt,U,v_old,resolvent,S)

	v = resolvent(mu + dt*S(U,dx) + v_old);

end
