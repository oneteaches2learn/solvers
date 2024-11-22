function my_pplane(x0,y0,Tend)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % M. Peszynska for MTH 4/581
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % plot phase portrait for x'=Ax and a trajectory with some I.C.
    % example (MUST uncomment one of lines [tode,yode] below
    % my_pplane(1,2,5)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [xx,yy]=meshgrid(-1:.1:1,-1:.1:1);

    % calculate slope fields when A = [-2 0; 0 -1]
    % define the functions
    dx = @(x,y)(-x); dy = @(x,y)(-2*y);

    % now actually calculate these on a grid
    dxx = dx(xx,yy); dyy = dy(xx,yy);

    % plot these slope fields
    quiver(xx,yy,dxx,dyy);

    % now find a trajectory
    tdum = 0;
    pause

    % EXPLORE two different ways to set up ode45 for t in [0,Tend]
    %UNCOMMENT one of the two lines
    %[tode,yode] = ode45(@linear_direct,[0,Tend],[x0;y0]);
    [tode,yode] = ode45(@linear_with_matrix,[0,Tend],[x0;y0]);

    hold on;
    plot(yode(:,1),yode(:,2));
    pause
    hold off;

    % plot in time
    plot(tode,yode(:,1),'k-',tode,yode(:,2),'r+');
    legend('y1','y2');
end
    
function dy = linear_direct(t,y)
    dy = zeros(2,1); %% this line is not mandatory
    dy = [-2*y(1,:);-y(2,:)];
end

function dy = linear_with_matrix(t,y)
    dy = zeros(2,1); %% this line is not mandatory
    A = [-2,0;0,-1];
    dy = A*y;
end
