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
            
            % Post-Processing
            %self = self.postProcess();
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
            % Assembles the global load vector F
            
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
        
        %% POST-PROCESSING
        function self = postProcess(self)
            % Performs post-processing tasks:
            % - Computes analytical solution
            % - Calculates L2 error
            % - Generates plots
            
            x = self.domain.x;
            u_numeric = self.solution;
            
            % Analytical solution (modify if using MMS)
            u_exact = (1 / pi^2) * sin(pi * x);
            
            % Store analytical solution
            self.post.u_exact = u_exact;
            
            % Calculate Absolute Error
            error = abs(u_numeric - u_exact);
            self.post.error = error;
            
            % Calculate L2 Error
            h = self.domain.h;
            L2_error = sqrt(sum(error.^2) * h);
            self.post.L2_error = L2_error;
            
            % Display L2 Error
            fprintf('L2 Error: %.5e\n', L2_error);
            
            % Plot Numerical and Analytical Solutions
            figure;
            plot(x, u_numeric, 'bo-', 'LineWidth', 1.5, 'MarkerFaceColor', 'b');
            hold on;
            plot(x, u_exact, 'r--', 'LineWidth', 2);
            xlabel('Spatial Coordinate x');
            ylabel('Solution u(x)');
            title('Finite Element Solution vs. Analytical Solution');
            legend('Numerical FEM Solution', 'Analytical Solution');
            grid on;
            hold off;
            
            % Plot Absolute Error
            figure;
            plot(x, error, 'k-o', 'LineWidth', 1.5, 'MarkerFaceColor', 'k');
            xlabel('Spatial Coordinate x');
            ylabel('|u_{numerical} - u_{exact}|');
            title('Absolute Error Between Numerical and Analytical Solutions');
            grid on;
        end
    end
end