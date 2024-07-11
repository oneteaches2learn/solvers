function [dt,nt,tGrid,init] = constructTimeStepping(time)

    dt    = time{1};
    nt    = time{2};
    init  = time{3};
    tGrid = linspace(0,nt*dt,nt+1);

end

