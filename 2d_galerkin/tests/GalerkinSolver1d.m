% GalerkinSolver1d.m
% A MATLAB class for performing 1D Galerkin Finite Element Analysis

classdef GalerkinSolver1d
    
    properties
        domain          % Structure containing domain parameters
        coefficients    % Structure containing physical coefficients
        solution        % Full solution vector including boundary nodes
    end
    
    properties (Hidden)
        tensors         % Hidden properties for internal matrix storage
        vectors         % Hidden properties for internal vector storage
    end
    
    methods
        %% CONSTRUCTOR
        function self = GalerkinSolver1d(domain, coefficients)
            % Constructor for GalerkinSolver1d
            %
            % Inputs:
            %   - domain: Structure with fields:
            %       M: Number of finite elements
            %       L: Length of the domain
            %       x: Grid points (column vector)
            %       h: Element size
            %       N: Number of interior nodes
            %   - coefficients: Structure with fields:
            %       k: Diffusion coefficient
            %       f: Function handle for source term f(x)
            %       r: Reaction coefficient
            
            if nargin == 2
                self.domain = domain;
                self.coefficients = coefficients;
            else
                error('Constructor requires domain and coefficients structures as inputs.');
            end
            
            % Assemble Stiffness Matrix and Load Vector
            self = self.assembleStiffnessMatrix();
            self = self.assembleMassMatrix(self.coefficients.r);
            self = self.assembleLoadVector();
            
            % Apply Boundary Conditions
            self = self.applyBoundaryConditions();
            
            % Solve the Linear System
            self = self.solveSystem();
        end
        
        %% ASSEMBLE STIFFNESS MATRIX
        function self = assembleStiffnessMatrix(self)
            % Assembles the global stiffness matrix K
            
            M = self.domain.M;
            N = self.domain.N;
            h = self.domain.h;
            k = self.coefficients.k;
            
            % Initialize global stiffness matrix
            K = sparse(N, N);
            
            % Local stiffness matrix for each element
            K_local = (k / h) * [1, -1; -1, 1];
            
            % Assemble global stiffness matrix
            for e = 1:M
                % Nodes associated with the current element
                n1 = e;     % Global node number (1 to M)
                n2 = e + 1; % Global node number (2 to M+1)
                
                % Map local nodes to global unknowns
                global_nodes = [n1, n2];
                global_unknowns = global_nodes - 1;  % Interior nodes are 2 to M
                
                % Assemble the stiffness matrix
                for i = 1:2
                    for j = 1:2
                        if (global_unknowns(i) >= 1) && (global_unknowns(i) <= N) && ...
                           (global_unknowns(j) >= 1) && (global_unknowns(j) <= N)
                            K(global_unknowns(i), global_unknowns(j)) = ...
                                K(global_unknowns(i), global_unknowns(j)) + K_local(i, j);
                        end
                    end
                end
            end
            
            % Store the assembled stiffness matrix
            self.tensors.K = K;
        end
        
        %% ASSEMBLE MASS MATRIX
        function self = assembleMassMatrix(self, c)
            % Assembles the global mass matrix R scaled by coefficient c
            
            M = self.domain.M;
            N = self.domain.N;
            h = self.domain.h;
            
            % Initialize global mass matrix
            R = sparse(N, N);
            
            % Local mass matrix for each element
            R_local = c * (h / 6) * [2, 1; 1, 2];
            
            % Assemble global mass matrix
            for e = 1:M
                % Nodes associated with the current element
                n1 = e;
                n2 = e + 1;
                
                % Map local nodes to global unknowns
                global_nodes = [n1, n2];
                global_unknowns = global_nodes - 1;  % Interior nodes
                
                % Assemble the mass matrix
                for i = 1:2
                    for j = 1:2
                        if (global_unknowns(i) >= 1) && (global_unknowns(i) <= N) && ...
                           (global_unknowns(j) >= 1) && (global_unknowns(j) <= N)
                            R(global_unknowns(i), global_unknowns(j)) = ...
                                R(global_unknowns(i), global_unknowns(j)) + R_local(i, j);
                        end
                    end
                end
            end
            
            % Store the assembled mass matrix
            self.tensors.R = R;
        end
        
        %% ASSEMBLE LOAD VECTOR
        function self = assembleLoadVector(self)
            % Assembles the global load vector F using the Trapezoidal Rule
            
            M = self.domain.M;
            N = self.domain.N;
            h = self.domain.h;
            x = self.domain.x;
            f = self.coefficients.f;
            
            % Initialize global load vector
            F = zeros(N, 1);
            
            % Local load vector for each element using Trapezoidal Rule
            for e = 1:M
                % Nodes associated with the current element
                n1 = e;     % Global node number (1 to M)
                n2 = e + 1; % Global node number (2 to M+1)
                
                % Compute local load vector
                F_local = (h / 6) * [2*f(x(n1)) + f(x(n2)); f(x(n1)) + 2*f(x(n2))];
                
                % Map local nodes to global unknowns
                global_nodes = [n1, n2];
                global_unknowns = global_nodes - 1;  % Interior nodes are 2 to M
                
                % Assemble the load vector
                for i = 1:2
                    if (global_unknowns(i) >= 1) && (global_unknowns(i) <= N)
                        F(global_unknowns(i)) = F(global_unknowns(i)) + F_local(i);
                    end
                end
            end
            
            % Store the assembled load vector
            self.vectors.F = F;
        end
        
        %% APPLY BOUNDARY CONDITIONS
        function self = applyBoundaryConditions(self)
            % Applies homogeneous Dirichlet boundary conditions
            % Placeholder for future extensions
        end
        
        %% SOLVE THE LINEAR SYSTEM
        function self = solveSystem(self)
            % Solves the linear system (K + R) * u_interior = F
            
            K = self.tensors.K;
            R = self.tensors.R;
            F = self.vectors.F;
            
            % Combine stiffness and mass matrices
            K_total = K + R;
            
            % Solve for interior node values
            u_interior = K_total \ F;
            
            % Construct the full solution including boundary nodes
            self.solution = [0; u_interior; 0];
        end
        
        %% EVALUATE FUNCTION ON GRID
        function uTrue_vector = evaluateFunction(self, uTrue)
            % Evaluates the input function on the grid in the domain.
            % Input can be a vector, function_handle, or symbolic function.
            %
            % Inputs:
            %   - uTrue: vector, function_handle, or symbolic function
            %
            % Outputs:
            %   - uTrue_vector: vector evaluated on the grid
            
            if isvector(uTrue) && isnumeric(uTrue)
                if length(uTrue) == length(self.domain.x)
                    uTrue_vector = uTrue(:); % Ensure column vector
                else
                    error('The length of the input vector does not match the number of nodes in the domain.');
                end

            elseif isa(uTrue, 'function_handle')
                uTrue_vector = uTrue(self.domain.x);
                if ~isvector(uTrue_vector) || length(uTrue_vector) ~= length(self.domain.x)
                    error('The function_handle provided does not return a vector matching the number of nodes in the domain.');
                end
                uTrue_vector = uTrue_vector(:); % Ensure column vector

            elseif isa(uTrue, 'sym') || isa(uTrue,'symfun')
                % Convert symbolic function to matlab function
                try
                    syms x
                    uTrue_func = matlabFunction(uTrue, 'Vars', x);
                    uTrue_vector = uTrue_func(self.domain.x);
                    if ~isvector(uTrue_vector) || length(uTrue_vector) ~= length(self.domain.x)
                        error('The symbolic function provided does not evaluate to a vector matching the number of nodes in the domain.');
                    end
                    uTrue_vector = uTrue_vector(:); % Ensure column vector
                catch
                    error('Failed to convert symbolic function to a MATLAB function handle.');
                end

            else
                error('uTrue must be a vector, function_handle, or symbolic function.');
            end
        end
        
        %% PLOT NUMERICAL AND ANALYTICAL SOLUTION
        function plot(self, uTrue)
            % Plots the numerical FEM solution.
            % If uTrue is provided, plots both numerical and true solutions.
            
            x = self.domain.x;
            u_numeric = self.solution;
            
            figure;
            hold on;
            
            % Determine Marker Indices
            num_nodes = length(x);
            if num_nodes > 40
                num_markers = 21;
                indices = round(linspace(1, num_nodes, num_markers));
            else
                indices = 1:num_nodes;
            end
            
            % Plot Numerical Solution with Markers
            plot(x, u_numeric, 'b-', 'LineWidth', 1.5, 'Marker', 'o', ...
                 'MarkerIndices', indices, 'MarkerFaceColor', 'b');
            
            if nargin > 1 && ~isempty(uTrue)
                uTrue_vector = self.evaluateFunction(uTrue);
                plot(x, uTrue_vector, 'r--', 'LineWidth', 1.5, 'Marker', 'o', ...
                     'MarkerIndices', indices, 'MarkerFaceColor', 'r');
                legend('Numerical FEM Solution', 'True Solution');
            else
                legend('Numerical FEM Solution');
            end
            
            xlabel('Spatial Coordinate x');
            ylabel('Solution u(x)');
            title('Finite Element Solution');
            grid on;
            hold off;
        end
        
        %% PLOT ABSOLUTE ERROR
        function plotError(self, uTrue)
            % Plots the absolute error between numerical and true solutions.
            
            if nargin < 2 || isempty(uTrue)
                error('plotError requires uTrue as an input.');
            end
            
            uTrue_vector = self.evaluateFunction(uTrue);
            x = self.domain.x;
            u_numeric = self.solution;
            error = abs(u_numeric - uTrue_vector);
            
            figure;
            hold on;
            
            % Determine Marker Indices
            num_nodes = length(x);
            if num_nodes > 40
                num_markers = 21;
                indices = round(linspace(1, num_nodes, num_markers));
            else
                indices = 1:num_nodes;
            end
            
            % Plot Absolute Error with Markers
            plot(x, error, 'k-', 'LineWidth', 1.5, 'Marker', 'o', ...
                 'MarkerIndices', indices, 'MarkerFaceColor', 'k');
            
            xlabel('Spatial Coordinate x');
            ylabel('|u_{numerical} - u_{True}|');
            title('Absolute Error Between Numerical and True Solutions');
            grid on;
            hold off;
        end
        
        %% COMPUTE L2 ERROR
        function L2_error = L2error(self, uTrue)
            % Computes the L2 norm of the error between numerical and true solutions.
            
            if nargin < 2 || isempty(uTrue)
                error('L2error requires uTrue as an input.');
            end
            
            uTrue_vector = self.evaluateFunction(uTrue);
            u_numeric = self.solution;
            error_vector = u_numeric - uTrue_vector;
            h = self.domain.h;
            L2_error = sqrt(sum(error_vector.^2) * h);
        end
    end
end