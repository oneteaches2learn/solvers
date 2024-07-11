clear all
% This code illustrates the use of change of coordinates to compute an integral
% over some triangle in a triangulated mesh. Triangle A, the mesh triangle, is
% defined by storing its coordinates as A.coord. These coordinates are used to
% create a square (3x3) matrix A.extended. The extended matrix is used to
% compute and store the area of A as A.area. A similar object, Ahat, is
% created, representin the reference triangle. The extended matrices of A and
% Ahat are used to define T.extended, from which T.transform and T.translate
% are stored. These entities represent the transformation:
%
% 		[x,y]' = T.transform * [xhat,yhat]' + T.translate
%
% where [x,y]' is a point in A and [xhat,yhat]' is a point in Ahat. The
% individual coordinate functions, stored as function_handles, that map x =
% x(xhat,yhat) and y = y(xhat,yhat) are also stored as T.x and T.y. Finally,
% T.J stores the jacobian determinant of T. 
%
% From all of the above, the integral of function f over A is computed in three
% ways: First, directly, by decomposing triangle A into two subregions and
% integrating using Matlab's integrate2 function. Second, by executing the
% change of coordinates using transform T and integrating over Ahat, again
% using Matlab's integrate2 function. Third, by again usin the change of
% coordinates, but this time integrating 'by hand'. What I mean is that I did
% the change of variables and integrated the resulting polynomial on paper. But
% then I am using Matlab to plug in the actual coefficients / numbers and do
% the final arithmetic. 
%
% The end result is that the three versions of the integral agree, as you would
% expect. Along the way, the numbers were used to generate an example in the
% documentation called BoundaryConditions.m
%
% author: Tyler Fara				Date: June 23, 2024
%-----------------------------------------------------------------------------%


% INPUTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% actual triangle
A.coord = [0.5 1; 0.644861 0.644861; 1 1];

% reference triangle
Ahat.coord = [1 0; 0 1; 0 0];

% desired function
f = @(x,y)(x.*(x-1).*y.*(y-1));


% TRANSFORMATION MATRIX %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% extend matrices
A.extended = [A.coord'; 1 1 1];
Ahat.extended = [Ahat.coord'; 1 1 1];

% triangle areas
A.area = abs(0.5 * det(A.extended));
Ahat.area = abs(0.5 * det(Ahat.extended));

% compute extended transformation matrix
T.extended = A.extended * Ahat.extended^-1;

% store rotation/scaling matrix and translation vector
T.transform = T.extended(1:2,1:2);
T.translate = T.extended(1:2,3);

% store x and y coordinate functions
T.x = @(xhat,yhat)(T.transform(1,1) * xhat + T.transform(1,2) * yhat + T.translate(1));
T.y = @(xhat,yhat)(T.transform(2,1) * xhat + T.transform(2,2) * yhat + T.translate(2));


% compute Jacobian determinant
T.J = det(T.transform(1:2,1:2));
Ahat.area / A.area;


% INTEGRALS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
q_actual = actual_integral(A.coord,f)
q_transf = abs(T.J) * transformed_integral(A.coord,f,T)


% PLOTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

plot_tri(A.coord,1)
plot_tri(Ahat.coord,2)


% 'by hand' calculation of integral
a = -0.5;
b = -0.3551;
c1 = (1/2) * a^2 * b^1 * (1/30);
c2 = (1/5) * a^0 * b^4 * (1/6);
c3 = (1/2) * a^0 * b^3 * (1/5);
c4 = (1/2) * a^1 * b^1 * (1/12);
c5 = (1/3) * a^0 * b^2 * (1/4);
c6 = (1/3) * a^2 * b^2 * (1/60);
c7 = (1/2) * a^1 * b^3 * (1/30);
c8 = (1/1) * a^1 * b^2 * (1/20);
c = (c1 + c2 + c3 + c4 + c5 + c6 + c7 + c8) * T.J;



% AUXILLIARY FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plot_tri(A,option)

	if option == 1
		color = "#0072BD";
		labels = {'aHat','bHat','cHat'};
	elseif option == 2
		color = "#7E2F8E";
		labels = {'a','b','c'};
	end

	A = [A; A(1,:)];

	hold on
	% plot edges
	for i = 1:3
		plot(A(i:i+1,1),A(i:i+1,2),'Color',color,'LineWidth',4);
	end

	% plot vertices
	plot(A(:,1),A(:,2),'.','Color',color,'MarkerSize',80);
	for i = 1:3
		text(A(i,1)+0.03,A(i,2)-0.04,labels{i},'FontSize',40);
	end
	hold off

end

function q = actual_integral(A,f)

	% sort coordinate matrix by x coordinate
	[Bx,I] = sort(A(:,1));
	By = A(:,2);
	By = By(I);
	A = [Bx,By];
	
	% store x-coordinates
	x1 = A(1,1);
	x2 = A(2,1);
	x3 = A(3,1);

	% compute function_handles for lower domain bounds
	y1 = line_computer(A(1,:),A(2,:));
	y2 = line_computer(A(2,:),A(3,:));

	% compute integrals
	q1 = integral2(f,x1,x2,y1,1);
	q2 = integral2(f,x2,x3,y2,1);
	
	q = q1 + q2;

end


function q = transformed_integral(A,f,T)
	
	f_T = @(xhat,yhat)(f(T.x(xhat,yhat),T.y(xhat,yhat)))
	y = line_computer([0,1],[1,0]);
	q = integral2(f_T,0,1,0,y);

end



function f = line_computer(a,b)

	coef = polyfit([a(1),b(1)],[a(2),b(2)],1);
	f = @(x)(coef(1) * x + coef(2));

end
