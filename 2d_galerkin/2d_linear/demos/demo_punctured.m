% DEMO_PUNCTURED Demonstration script for punctured 2D domains
clear all

%{
% [DEMO 1] Create and plot a simple domain
% domain bounds
xLim_dom = [0 2];
yLim_dom = [0 1];

% create Domain2d object
dom = GalerkinAssembler2d_poisson.assembleDomainGeometry(xLim_dom,yLim_dom);

% plot the domain
dom.boundary.plot;
pause();
close();


% [DEMO 2] Create and plot a punctured domain
% domain bounds
xLim_dom = [0 2];
yLim_dom = [0 1];

% reference cell bounds
xLim_Y = [0 1];
yLim_Y = [0 1];

% inclusion parameters
epsilon = 1;
incRatio = pi / 2;

% create Inclusion2d object
inc_circle = Inclusion2d_circle(xLim_Y,yLim_Y,incRatio);
inc_square = Inclusion2d_square(xLim_Y,yLim_Y,incRatio);

% create Domain2d_punctured object
dom_circle = GalerkinAssembler2d_poisson.assembleDomainGeometry(xLim_dom,yLim_dom,inc_circle,epsilon);
dom_square = GalerkinAssembler2d_poisson.assembleDomainGeometry(xLim_dom,yLim_dom,inc_square,epsilon);

% plot the domains
figure;
subplot(2,1,1);
dom_circle.boundary.plot;
subplot(2,1,2);
dom_square.boundary.plot;
pause();
close();

% plot the domain, with labels
figure;
subplot(2,1,1);
dom_circle.boundary.plot("EdgeLabels","on","FaceLabels","on");
subplot(2,1,2);
dom_square.boundary.plot("EdgeLabels","on","FaceLabels","on");
pause();
close();


% [DEMO 3] Create and plot a punctured domain with smaller inclusions
% domain bounds
xLim_dom = [0 2];
yLim_dom = [0 1];

% reference cell bounds
xLim_Y = [0 1];
yLim_Y = [0 1];

% inclusion parameters
epsilon = 1;
incRatio = pi / 4; % <~~~ incRatio = |delta Q| / |Y|

% create Inclusion2d object
inc_circle = Inclusion2d_circle(xLim_Y,yLim_Y,incRatio);
inc_square = Inclusion2d_square(xLim_Y,yLim_Y,incRatio);

% create Domain2d_punctured object
dom_circle = GalerkinAssembler2d_poisson.assembleDomainGeometry(xLim_dom,yLim_dom,inc_circle,epsilon);
dom_square = GalerkinAssembler2d_poisson.assembleDomainGeometry(xLim_dom,yLim_dom,inc_square,epsilon);

% plot the domains
figure;
subplot(2,1,1);
dom_circle.plotGeometry;
subplot(2,1,2);
dom_square.plotGeometry;
pause();
close();


% [DEMO 4] Create and plot a punctured domain with 200 inclusions
% domain bounds
xLim_dom = [0 2];
yLim_dom = [0 1];

% reference cell bounds
xLim_Y = [0 1];
yLim_Y = [0 1];

% inclusion parameters
epsilon = 0.1;
incRatio = pi / 2; % <~~~ incRatio = |delta Q| /

% create Inclusion2d object
inc_circle = Inclusion2d_circle(xLim_Y,yLim_Y,incRatio);
inc_square = Inclusion2d_square(xLim_Y,yLim_Y,incRatio);

% create Domain2d_punctured object
dom_circle = GalerkinAssembler2d_poisson.assembleDomainGeometry(xLim_dom,yLim_dom,inc_circle,epsilon);
dom_square = GalerkinAssembler2d_poisson.assembleDomainGeometry(xLim_dom,yLim_dom,inc_square,epsilon);

% plot the domains
figure;
subplot(2,1,1);
dom_circle.plotGeometry;
subplot(2,1,2);
dom_square.plotGeometry;
pause();
close();


% [DEMO 5a] Mesh a punctured domain
% domain bounds
xLim_dom = [0 2];
yLim_dom = [0 1];

% reference cell bounds
xLim_Y = [0 1];
yLim_Y = [0 1];

% inclusion parameters
epsilon = 0.5;
incRatio = pi / 2; % <~~~ incRatio = |delta Q| / |Y|

% mesh parameters 
base = 2;
p_coarse = 2;
p_fine = 5;

% create Inclusion2d object
inc_circle = Inclusion2d_circle(xLim_Y,yLim_Y,incRatio);
inc_square = Inclusion2d_square(xLim_Y,yLim_Y,incRatio);

% create Domain2d_punctured object
dom_circle = GalerkinAssembler2d_poisson.assembleDomainGeometry(xLim_dom,yLim_dom,inc_circle,epsilon);
dom_square = GalerkinAssembler2d_poisson.assembleDomainGeometry(xLim_dom,yLim_dom,inc_square,epsilon);

% mesh the domains
dom_circle_coarse = GalerkinAssembler2d_poisson.assembleMesh(dom_circle,p_coarse,base);
dom_circle_fine = GalerkinAssembler2d_poisson.assembleMesh(dom_circle,p_fine,base);
dom_square_coarse = GalerkinAssembler2d_poisson.assembleMesh(dom_square,p_coarse,base);
dom_square_fine = GalerkinAssembler2d_poisson.assembleMesh(dom_square,p_fine,base);

% plot the mesh
figure;
subplot(2,2,1);
dom_circle_coarse.mesh.plot;
title('Coarse Mesh, Circular Inclusions');
subplot(2,2,2);
dom_square_coarse.mesh.plot;
title('Coarse Mesh, Square Inclusions');
subplot(2,2,3);
dom_circle_fine.mesh.plot;
title('Fine Mesh, Circular Inclusions');
subplot(2,2,4);
dom_square_fine.mesh.plot;
title('Fine Mesh, Square Inclusions');
set(gcf,'Position',[100 100 1200 800]);
pause();
close();


% [DEMO 6] Demo mesh properties
% domain bounds
xLim_dom = [0 1];
yLim_dom = [0 1];

% reference cell bounds
xLim_Y = [0 1];
yLim_Y = [0 1];

% inclusion parameters
epsilon = 1;
incRatio = pi / 2; % <~~~ incRatio = |delta Q| / |Y|

% mesh parameters 
p = 1;
base = 2;

% create and mesh domain
inc = Inclusion2d_circle(xLim_Y,yLim_Y,incRatio);
dom = GalerkinAssembler2d_poisson.assembleDomainGeometry(xLim_dom,yLim_dom,inc,epsilon);
dom = GalerkinAssembler2d_poisson.assembleMesh(dom,p,base);

% display mesh properties
mesh = dom.mesh
pause();

% display mesh nodes
nodes = dom.mesh.nodes

% plot mesh nodes on domain
figure;
dom.mesh.plot;
hold on;
plot(dom.mesh.nodes(:,1),dom.mesh.nodes(:,2),'k.','MarkerSize',48);
title('Mesh Nodes');
hold off;
pause();
close();

% display mesh midpoints
midpoints = dom.mesh.midpoints

% plot mesh midpoints on domain
figure;
dom.mesh.plot;
hold on;
plot(dom.mesh.midpoints(:,1),dom.mesh.midpoints(:,2),'k.','MarkerSize',48);
title('Mesh Midpoints');
hold off;
pause();
close();

% display mesh centroids
centroids = dom.mesh.centroids

% plot mesh centroids on domain
figure;
dom.mesh.plot;
hold on;
plot(dom.mesh.centroids(:,1),dom.mesh.centroids(:,2),'k.','MarkerSize',48);
title('Mesh Centroids');
hold off;
pause();
close();

% display mesh node numbers
dom.mesh.plot("NodeLabels","on");
title('Mesh Nodes Labels');
pause();
close();

% display mesh element numbers
dom.mesh.plot("ElementLabels","on");
title('Mesh Element Labels');
pause();
close();

% display mesh areas
areas = dom.mesh.areas

% plot 
figure;
dom.mesh.plot;
hold on;
text(dom.mesh.centroids(:,1),dom.mesh.centroids(:,2),string(dom.mesh.areas),'FontSize',12,'Color','k','HorizontalAlignment','center');
hold off;
title('Element Areas');
pause();
close();

% display mesh edge lengths
edgeLengths = dom.mesh.edgeLengths

% plot mesh edge lengths on domain
figure;
dom.mesh.plot;
hold on;
text(dom.mesh.midpoints(:,1),dom.mesh.midpoints(:,2),string(round(dom.mesh.edgeLengths,2)),'FontSize',12,'Color','k','HorizontalAlignment','center');
hold off;
title('Edge Lengths');
pause();
close();

% display mesh quality
quality = dom.mesh.quality
quality_byElement = dom.mesh.quality_byElement

% plot mesh quality on domain
figure;
dom.mesh.plot;
hold on;
text(dom.mesh.centroids(:,1),dom.mesh.centroids(:,2),string(round(dom.mesh.quality_byElement,2)),'FontSize',12,'Color','k','HorizontalAlignment','center');
hold off;
title('Element Quality');
pause();
close();


%{
% I think delete this one?
% [DEMO 7] Plot a function on a domain
% domain bounds
xLim_dom = [0 2];
yLim_dom = [0 1];

% reference cell bounds
xLim_Y = [0 1];
yLim_Y = [0 1];

% inclusion parameters
epsilon = 0.5;
incRatio = pi / 2; % <~~~ incRatio = |delta Q| / |Y|

% mesh parameters 
p = 4;
base = 2;

% create Inclusion2d object
inc_circle = Inclusion2d_circle(xLim_Y,yLim_Y,incRatio);
inc_square = Inclusion2d_square(xLim_Y,yLim_Y,incRatio);

% create Domain2d_punctured object
dom_circle = GalerkinAssembler2d_poisson.assembleDomainGeometry(xLim_dom,yLim_dom,inc_circle,epsilon);
dom_square = GalerkinAssembler2d_poisson.assembleDomainGeometry(xLim_dom,yLim_dom,inc_square,epsilon);

% mesh the domains
dom_circle = GalerkinAssembler2d_poisson.assembleMesh(dom_circle,p,base);
dom_square = GalerkinAssembler2d_poisson.assembleMesh(dom_square,p,base);

% plot the mesh
figure;
subplot(2,1,1);
dom_circle.mesh.plot;
subplot(2,1,2);
dom_square.mesh.plot;
pause();
close();


% I think delete this one?
% [DEMO 5b] Mesh a punctured domain and plot with labels
% domain bounds
xLim_dom = [0 2];
yLim_dom = [0 1];

% reference cell bounds
xLim_Y = [0 1];
yLim_Y = [0 1];

% inclusion parameters
epsilon = 1;
incRatio = pi / 2; % <~~~ incRatio = |delta Q| / |Y|

% mesh parameters 
p = 2;
base = 2;

% create Inclusion2d object
inc_circle = Inclusion2d_circle(xLim_Y,yLim_Y,incRatio);
inc_square = Inclusion2d_square(xLim_Y,yLim_Y,incRatio);

% create Domain2d_punctured object
dom_circle = GalerkinAssembler2d_poisson.assembleDomainGeometry(xLim_dom,yLim_dom,inc_circle,epsilon);
dom_square = GalerkinAssembler2d_poisson.assembleDomainGeometry(xLim_dom,yLim_dom,inc_square,epsilon);

% mesh the domains
dom_circle = GalerkinAssembler2d_poisson.assembleMesh(dom_circle,p,base);
dom_square = GalerkinAssembler2d_poisson.assembleMesh(dom_square,p,base);

% plot the mesh, with labels
figure;
subplot(2,1,1);
dom_circle.mesh.plot("NodeLabels","on","ElementLabels","on");
subplot(2,1,2);
dom_square.mesh.plot("NodeLabels","on","ElementLabels","on");
pause();
close(); 
%}


% [DEMO 6] Evaluate a function_handle on a domain
% domain bounds
xLim_dom = [0 1];
yLim_dom = [0 1];

% reference cell bounds
xLim_Y = [0 1];
yLim_Y = [0 1];

% inclusion parameters
epsilon = 1;
incRatio = pi / 2; % <~~~ incRatio = |delta Q| / |Y|

% mesh parameters 
base = 2;
p = 1;

% create and mesh a domain with inclusions
inc = Inclusion2d_circle(xLim_Y,yLim_Y,incRatio);
dom = GalerkinAssembler2d_poisson.assembleDomainGeometry(xLim_dom,yLim_dom,inc,epsilon);
dom = GalerkinAssembler2d_poisson.assembleMesh(dom,p,base);

% define a function_handle to plot
f_han = @(x,y) sin(-pi*x) .* cos(pi*y);

% evalauate the function on the mesh nodes
f_nod = dom.function2nodes(f_han)
pause();
f_mid = dom.function2midpoints(f_han)
pause();
f_cen = dom.function2centroids(f_han)
pause();


% [DEMO 7] Evaluate a function_handle on a domain
% domain bounds
xLim_dom = [0 2];
yLim_dom = [0 1];

% reference cell bounds
xLim_Y = [0 1];
yLim_Y = [0 1];

% inclusion parameters
epsilon = 0.5;
incRatio = pi / 2; % <~~~ incRatio = |delta Q| / |Y|

% mesh parameters 
base = 2;
p = 5;

% create and mesh a domain with inclusions
inc = Inclusion2d_circle(xLim_Y,yLim_Y,incRatio);
dom = GalerkinAssembler2d_poisson.assembleDomainGeometry(xLim_dom,yLim_dom,inc,epsilon);
dom = GalerkinAssembler2d_poisson.assembleMesh(dom,p,base);

% define a function_handle to plot
f_han = @(x,y) sin(-pi*x) .* cos(pi*y);

% evalauate the function on the mesh nodes
f_nod = dom.function2nodes(f_han);
f_mid = dom.function2midpoints(f_han);
f_cen = dom.function2centroids(f_han);


% 3d plot of the function handle
dom.plot_function(f_han);
title('3d Plot using Function Handle');
pause();
close();

% 2d patch plot of the function handle
dom.plotPatch_function(f_han);
title('2d Plot using Function Handle');
pause();
close();

% 2d plot of nodal vector
dom.plotPatch_function(f_nod);
title('2d Plot using Nodal Vector');
pause();
close();

% try to plot with midpoint vector
try dom.plotPatch_function(f_mid);
catch ME
    disp('Error plotting with midpoint vector:');
    disp(ME.message);
    pause();
end
%}


%{ 
% this one's not quite ready as I still need to check everything and then make
% it output the results as tables 
% [DEMO 9] Quadrature on a punctured domain
% domain bounds
xLim_dom = [0 1];
yLim_dom = [0 1];

% reference cell bounds
xLim_Y = [0 1];
yLim_Y = [0 1];

% inclusion parameters
epsilon = 1;
incRatio = 4/3; % <~~~ incRatio = |delta Q| / |Y|

% mesh parameters
base = 2;
p = {};
for i = 1:6, p{i} = i; end

% create and mesh domains with successfively finer meshes
inc = Inclusion2d_square(xLim_Y,yLim_Y,incRatio);
dom_geo = GalerkinAssembler2d_poisson.assembleDomainGeometry(xLim_dom,yLim_dom,inc,epsilon);
dom = {};
for i = 1:6
    tic;
    dom{i} = GalerkinAssembler2d_poisson.assembleMesh(dom_geo,p{i},base);
    meshTime{i} = toc;
end

% define function to integrate
f_han = @(x,y) sin(pi*x) .* sin(pi*y);

% plot function on finest mesh
dom{6}.plot_function(f_han);
title('Function to Integrate');

sol = 3 / pi^2;

% compute three point quadrature
ThreePt = {};
for i = 1:6
    % compute three point quadrature
    tic; ThreePt{i} = dom{i}.threePointQuadrature(f_han);
    threePtTime{i} = toc;

    % compute error
    err{i} = abs(sol - ThreePt{i});

    % compute rate
    if i == 1, rate{i} = NaN;
    else, rate{i-1} = log(err{i-1} / err{i}) / log(2); end
end

% compute nodal quadrature
Nodal = {};
for i = 1:6
    % compute nodal quadrature
    tic; Nodal{i} = dom{i}.nodalQuadrature(f_han);
    nodalTime{i} = toc;

    % compute error
    err{i} = abs(sol - Nodal{i});

    % compute rate
    if i == 1, rate{i} = NaN;
    else, rate{i-1} = log(err{i-1} / err{i}) / log(2); end
end

% compute centroid quadrature
Centroid = {};
for i = 1:6
    % compute centroid quadrature
    tic; Centroid{i} = dom{i}.centroidQuadrature(f_han);
    centroidTime{i} = toc;  

    % compute error
    err{i} = abs(sol - Centroid{i});

    % compute rate
    if i == 1, rate{i} = NaN;
    else, rate{i-1} = log(err{i-1} / err{i}) / log(2); end
end
%}

%{
% Y bounds
xLim_Y = [0 1];
yLim_Y = [0 1];

% number of inclusions
epsilon = 1;
incRatio = pi / 2; % <~~~ incRatio = |delta Q| / |Y|
%}


%{
% mms parameters
p = 2;
base = 2;

% specify coefficients
k = 2 + sin(x(1) + x(2));
r = 2 + sin(x(1) + x(2));

% specify source
f = 1;

% specify BCs
bTypes_outer = 'DDDD';
bTypes_inner = 'R';

% specify Dirichlet conditions
u_D = 0;

% specify Neumann conditions
u_N = 1;

% specify Robin conditions
beta = 1;
u_R = 1;



% MMS TEST %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Poisson Trial Begun\n')

% assemble domain
fprintf(' Contructing Domain:'), tic
	inc = Inclusion2d_circle(xLim_Y,yLim_Y,incRatio);
	%inc = Inclusion2d_square(xLim_Y,yLim_Y,incRatio);
	%dom = GalerkinAssembler2d_poisson.assembleDomainGeometry(xLim_dom,yLim_dom,inc,epsilon);
	%dom = GalerkinAssembler2d_poisson.assembleDomainGeometry(xLim_dom,yLim_dom);
	dom = Domain2d(xLim_dom,yLim_dom);
	%dom = dom.add_yline;
	%dom = dom.inclusionsOFF;
	dom = GalerkinAssembler2d_poisson.assembleBoundary(dom,bTypes_outer,u_D,u_N,beta,u_R,bTypes_inner); 
	dom = GalerkinAssembler2d_poisson.assembleMesh(dom,p,base); 
executionTime = toc; 
fprintf(' %f s\n',executionTime)

% assemble other parameters
fprintf(' Assembling Parameters:'), tic
	auxfun = GalerkinAssembler2d_poisson.assembleCoefficients(k,r,f);
executionTime = toc; 
fprintf(' %f s\n',executionTime)

% run trial
fprintf(' Generating Solution:'), tic
	prob = GalerkinSolver2d_poisson(dom,auxfun);
executionTime = toc; 
fprintf(' %f s\n',executionTime)
%}