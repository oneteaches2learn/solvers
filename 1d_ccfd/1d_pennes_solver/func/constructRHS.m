function Q = constructRHS(dx,xc,t,dt,f,C,theta,Por,nsol)
%CONSTRUCTRHS directs the construction of the right hand side.

    F = computeFunction(f,xc,t);
    Q = dt*F.*dx + dt*theta.*C.*dx + Por.*dx.*nsol;
    
end


