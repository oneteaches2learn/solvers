% GalerkinSolver1d.m
% A MATLAB class for performing 1D Galerkin Finite Element Analysis
% Author: [Your Name]
% Date: [Date]

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
            
            if nargin == 2
                self.domain = domain;
                self.coefficients = coefficients;
            else
                error('Constructor requires domain and coefficients structures as inputs.');
            end
            
            % Assemble Stiffness Matrix and Load Vector
            self = self.assembleStiffnessMatrix();
            self = self.assembleLoadVector();
            
            % Apply Boundary Conditions
            self = self.applyBoundaryConditions();
            
            % Solve the Linear System
            self = self.solveSystem();
            
            % Removed postProcess call
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
            %
            % For homogeneous Dirichlet BCs, and since boundary nodes are excluded
            % from the unknowns, no modifications to K and F are necessary.
            % If non-homogeneous BCs are required, modify F accordingly.
            
            % Placeholder for future extensions (e.g., non-homogeneous BCs)
        end
        
        %% SOLVE THE LINEAR SYSTEM
        function self = solveSystem(self)
            % Solves the linear system K * u_interior = F
            
            K = self.tensors.K;
            F = self.vectors.F;
            
            % Solve for interior node values
            u_interior = K \ F;
            
            % Construct the full solution including boundary nodes
            self.solution = [0; u_interior; 0];
        end
        
        %% PLOT NUMERICAL AND ANALYTICAL SOLUTION
        function plot(self, uTrue)
            % Plots the numerical FEM solution.
            % If uTrue is provided, plots both numerical and true solutions.
            %
            % Inputs:
            %   - uTrue (optional): Vector containing the true solution.
            
            x = self.domain.x;
            u_numeric = self.solution;
            
            figure;
            plot(x, u_numeric, 'bo-', 'LineWidth', 1.5, 'MarkerFaceColor', 'b');
            hold on;
            
            if nargin > 1 && ~isempty(uTrue)
                plot(x, uTrue, 'r--', 'LineWidth', 2);
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
            %
            % Inputs:
            %   - uTrue: Vector containing the true solution.
            %
            % Usage:
            %   fem_solver.plotError(uTrue);
            
            if nargin < 2 || isempty(uTrue)
                error('plotError requires uTrue as an input.');
            end
            
            x = self.domain.x;
            u_numeric = self.solution;
            error = abs(u_numeric - uTrue);
            
            figure;
            plot(x, error, 'k-o', 'LineWidth', 1.5, 'MarkerFaceColor', 'k');
            xlabel('Spatial Coordinate x');
            ylabel('|u_{numerical} - u_{True}|');
            title('Absolute Error Between Numerical and True Solutions');
            grid on;
        end
        
        %% COMPUTE L2 ERROR
        function L2_error = L2error(self, uTrue)
            % Computes the L2 norm of the error between numerical and true solutions.
            %
            % Inputs:
            %   - uTrue: Vector containing the true solution.
            %
            % Outputs:
            %   - L2_error: Scalar value representing the L2 norm of the error.
            %
            % Usage:
            %   L2_error = fem_solver.L2error(uTrue);
            
            if nargin < 2 || isempty(uTrue)
                error('L2error requires uTrue as an input.');
            end
            
            u_numeric = self.solution;
            error_vector = u_numeric - uTrue;
            h = self.domain.h;
            L2_error = sqrt(sum(error_vector.^2) * h);
        end
    end
end