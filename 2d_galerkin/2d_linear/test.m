% Define the problem domain and mesh
L = 1; % Length of the square domain
n = 10; % Number of elements along one dimension
h = L / n; % Element size
[x, y] = meshgrid(0:h:L, 0:h:L); % Generate grid points
nodes = [x(:), y(:)]; 
numNodes = size(nodes, 1);

% Define the connectivity matrix
elements = delaunay(nodes(:,1), nodes(:,2));
numElements = size(elements, 1);

% Initialize the stiffness matrix and load vector
K = sparse(numNodes, numNodes);
F = zeros(numNodes, 1);

% Define the source term (right-hand side of Poisson equation)
f = @(x, y) 1; 

% Assemble the stiffness matrix and load vector
for e = 1:numElements
    % Get the nodes of the current element
    elementNodes = elements(e, :);
    coords = nodes(elementNodes, :);
    
    % Compute the element stiffness matrix and load vector
    [Ke, Fe] = elementStiffness(coords, f);
    
    % Assemble into the global stiffness matrix and load vector
    K(elementNodes, elementNodes) = K(elementNodes, elementNodes) + Ke;
    F(elementNodes) = F(elementNodes) + Fe;
end

% Apply periodic boundary conditions
% Identify nodes on the left and right boundaries
leftNodes = find(nodes(:,1) == 0);
rightNodes = find(nodes(:,1) == L);

% Identify nodes on the bottom and top boundaries
bottomNodes = find(nodes(:,2) == 0);
topNodes = find(nodes(:,2) == L);

%{
% Map right boundary nodes to left boundary nodes
for i = 1:length(leftNodes)
    K(rightNodes(i), :) = K(rightNodes(i), :) + K(leftNodes(i), :);
    K(:, rightNodes(i)) = K(:, rightNodes(i)) + K(:, leftNodes(i));
    F(rightNodes(i)) = F(rightNodes(i)) + F(leftNodes(i));
end

% Map top boundary nodes to bottom boundary nodes
for i = 1:length(bottomNodes)
    K(topNodes(i), :) = K(topNodes(i), :) + K(bottomNodes(i), :);
    K(:, topNodes(i)) = K(:, topNodes(i)) + K(:, bottomNodes(i));
    F(topNodes(i)) = F(topNodes(i)) + F(bottomNodes(i));
end
%}

% Remove duplicate rows and columns
duplicates = union(leftNodes,bottomNodes);
K(duplicates,:) = [];
K(:,duplicates) = [];
F(duplicates) = [];

% Solve the linear system
U = K \ F;

% Reshape the solution for visualization
U_grid = reshape(U, n, n);

% Plot the solution
figure;
surf(x(1:end-1, 1:end-1), y(1:end-1, 1:end-1), U_grid);
title('Solution of the Poisson Problem with Periodic Boundary Conditions');
xlabel('x');
ylabel('y');
zlabel('u(x,y)');

% Function to compute the element stiffness matrix and load vector
function [Ke, Fe] = elementStiffness(coords, f)
    % Gauss quadrature points and weights for a triangle
    gaussPoints = [1/3, 1/3];
    gaussWeights = 1/2;
    
    % Shape functions and their gradients
    N = @(xi, eta) [1 - xi - eta, xi, eta];
    dN_dxi = [-1, 1, 0];
    dN_deta = [-1, 0, 1];
    
    % Jacobian matrix and its determinant
    J = [dN_dxi; dN_deta] * coords;
    detJ = det(J);
    
    % Gradient of shape functions in physical coordinates
    invJ = inv(J);
    dN_dx = invJ * [dN_dxi; dN_deta];
    
    % Element stiffness matrix
    Ke = (dN_dx' * dN_dx) * detJ * gaussWeights;
    
    % Element load vector
    Fe = N(gaussPoints(1), gaussPoints(2))' * f(mean(coords(:,1)), mean(coords(:,2))) * detJ * gaussWeights;
end