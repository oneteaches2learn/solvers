classdef Domain1d

    properties
        domainArea
        mesh
        boundary
        time
        xLim
    end

    methods
        % CONSTRUCTOR
        function self = Domain1d(x)

            if nargin == 1

                % create boundary object
                %self.boundary = Boundary1d(x);

                % store variables
                self.xLim = x;

                % compute area
                self.domainArea = x(2) - x(1);

            end
        end


        % SETTERS
        function self = set.mesh(self,mesh)
        % SET.MESH Set mesh property of the Domain1d object.
        % USAGE:
        %   self = set.mesh(self,mesh)
        % INPUTS:
        %   mesh - Mesh1d object

            % validate and assign input
            validateattributes(mesh, {'Mesh1d'}, {'scalar'});
            self.mesh = mesh;

            % if boundary is configured, assign nodes
            if ~isempty(self.boundary.bTypes) ...
                && length(self.boundary.bTypes) == self.boundary.nEdges

                self.boundary = self.boundary.mesh2nodes(mesh);

            else 
                warning(['Boundary1d object is not fully configured. ' ...
                    'Boundary nodes have not been assigned from the mesh.']);
            end

        end

        function self = set.time(self,varargin)
        % SET.TIME Set time property of the Domain1d object.
        % USAGE:
        %   self = set.time(self,time)
        %   self = set.time(self,T,dt)
        %   self = set.time(self,T,dt,eq)
        % INPUTS:
        %   time - Time1d object
        %   T    - total time
        %   dt   - time step size
        %   eq   - (optional) equilibrium status (true/false)

            % if input is a Time1d object
            if isa(varargin{1},'TimeStepping') 
                self.time = varargin{1};

            % else if inputs are time parameters
            elseif nargin >= 3
                T = varargin{1};
                dt = varargin{2};

                % if equilibrium status is provided
                if nargin == 4
                    eq = varargin{3};
                    self.time = Time1d(T,dt,eq);
                else
                    self.time = Time1d(T,dt);
                end
            end
        end


        % PLOTTERS
        function plot(self)

            % plot line from xLim(1) to xLim(2)
            plot([self.xLim(1), self.xLim(2)], [0, 0], 'k-', 'LineWidth', 2);
            hold on;

            % plot mesh nodes
            plot(self.mesh.nodes, zeros(size(self.mesh.nodes)), 'bo', 'MarkerSize', 8, 'MarkerFaceColor', 'b'); 

            % plot boundary points
            plot(self.xLim(1), 0, 'bs', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
            plot(self.xLim(2), 0, 'bs', 'MarkerSize', 10, 'MarkerFaceColor', 'r');

            % configure plot
            title('1D Domain Mesh');
            xlabel('x');
            ylim([-0.1,0.1]);
            yticks([]);
            grid on;

            % make the plot very narrow in the y direction
            axis tight;
            f = gcf;
            f.Position(4) = 150; 
            
            hold off;
        end

        function h = plot_function(self,func,timeStep)
        % PLOT_FUNCTION Plot a function over the domain.
        % INPUTS:
        %   func - function to be plotted; can be function handle, a symfun,
        %           a sym, a constant (i.e. a lone double), or a vector of doubles 
        %   timeStep - (optional) time step at which to evaluate the function
        % OUTPUTS:
        %   h    - handle to the plot 

            % STEP 1: time-step handling
            % Domain1d has time-stepping, but no time step is provided
            if nargin < 3 && ~isempty(self.time)
                tVal = self.time.tGrid(end);
            
            % Domain1d has time-stepping, and a time step is provided
            elseif nargin == 3 && ~isempty(self.time)
                tVal = self.time.tGrid(timeStep);

            % Domain1d has no time-stepping
            else
                timeStep = 1;
                tVal = 0;
            end 

            % if Domain1d has no time-stepping but a time step is provided, warn user
            if nargin == 3 && isempty(self.time)
                warning('Domain1d has no time-stepping, so the provided time step will be ignored.');
            end
        

            % STEP 2: function evaluation
            % store mesh nodes
            xVals = self.mesh.nodes;

            % if func is a vector of doubles, set it as yVals directly
            if isa(func,'double') && length(func) > 1
                if length(func) ~= self.mesh.nNodes
                    error('If func is a vector of doubles, its length must equal the number of mesh nodes.');
                end
                yVals = func;

            % else wrap and evaluate func
            else
                func = Domain1d.wrapFunction(func);
                yVals = func(xVals,tVal);

            end


            % STEP 3: plot and configure plot
            % create plot
            h = plot(xVals, yVals, 'LineWidth', 2);

            % configure plot
            title('Function over 1D Domain');
            xlabel('x');
            ylabel('f(x)');
            grid off;

            % if nNodes is small, plot evals on mesh nodes
            if self.mesh.nNodes <= 40
                hold on;

                % plot interior nodes in blue
                plot(self.mesh.nodes, yVals, 'bo', 'MarkerSize', ...
                    8, 'MarkerFaceColor', 'b');

                % plot boundary nodes in red
                plot(self.xLim(1), yVals(1), 'rs', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
                plot(self.xLim(2), yVals(end), 'rs', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
                hold off;
            end

        end

        function animate_function(self,func)
        % ANIMATE_FUNCTION Create an animation of a function over time.
        % INPUTS:
        %   func - function to be animated; can be a function handle, a symfun,
        %           a sym, a constant (i.e. a lone double), or a matrix of doubles
        %
        % NOTES:
        %   (1) If func is a matrix of doubles, its size must be nSpaceEvals x nTimeEvals.
        %
        %   (2) This method requires that the Domain1d object has a TimeStepping
        %   object assigned to its time property.
        %
        %   (3) The animation is created by evaluating the function to form an
        %   nSpaceEvals x nTimeEvals matrix. Then columns of this matrix are
        %   passed to the plot_function method to create the animation.
        %
        %   (4) Through experimentation, the mean time to plot a single time step
        %   is around 0.013 seconds. This value is used to estimate the total
        %   time of the animation. If the estimated total time is less than
        %   around 3 seconds, a pause buffer is added to slow down the animation
        %   so that it is visible to the user. If the estimated total time is
        %   more than around 3 seconds, the number of time steps plotted is
        %   reduced to around 200 to speed up the animation.
        
            
            % STEP 1: time-step handling
            % Domain1d has time-stepping
            if ~isempty(self.time)
                tVals = self.time.tGrid;
            
            % Domain1d has no time-stepping
            else
                error('Domain1d has no time-stepping. Animation cannot be created.');
            end 


            % STEP 2: function evaluation
            % store mesh nodes
            xVals = self.mesh.nodes;

            % if func is a vector of doubles, set it as yVals directly
            if isa(func,'double') && length(func) > 1
                if size(func,1) ~= self.mesh.nNodes || size(func,2) ~= self.time.N_t
                    error(['If func is a matrix of doubles, its size must be nSpaceEvals x nTimeEvals, ' ...
                        'where nSpaceEvals = %d and nTimeEvals = %d.'], self.mesh.nNodes, self.time.N_t);
                end
                yVals = func;

            % else wrap and evaluate func
            else
                func = Domain1d.wrapFunction(func);
                yVals = func(xVals,tVals);

            end


            % STEP 3: plot first time step
            % get data
            yMax = max(yVals,[],'all');
            yMin = min(yVals,[],'all');
            yLim = [yMin - 0.1*abs(yMin), yMax + 0.1*abs(yMax)];

            % generate plot
            h = self.plot_function(yVals(:,1));

            % configure plot
            titleStr = sprintf('Function over 1D Domain, t = %.2f', self.time.tGrid(1));
            title(titleStr);
            ylim(yLim);

            % pause to show first plot
            pause();


            % STEP 4: create animation
            meanStepTime = 0.013;
            estimatedTotalTime = meanStepTime * self.time.N_t;

            % if total time is less than ~3 seconds, add pause buffer
            if estimatedTotalTime <= 3
                pauseBuffer = (3 - estimatedTotalTime) / self.time.N_t;
                steps = [1:self.time.N_t];

            % if more than ~3 seconds, reduce number of steps to around 200
            else
                steps = [1:ceil(self.time.N_t / 200):self.time.N_t, self.time.N_t];
                pauseBuffer = 0;
            end

            % loop through time steps
            for n = steps

                % generate plot
                h = self.plot_function(yVals(:,n));

                % configure plot
                titleStr = sprintf('Function over 1D Domain, t = %.2f', self.time.tGrid(n));
                title(titleStr);
                ylim(yLim);

                % pause to create animation effect
                pause(pauseBuffer);
            end

        end


        % UTILITY FUNCTIONS
        function vals = function2nodes(self,func)
        % FUNCTION2NODES Evaluate a function at the mesh nodes.
        % INPUTS:
        %   func - function to be evaluated (can be function handle, symfun,
        %           sym, or constant i.e. double) 
        % OUTPUTS:
        %   vals - function values at mesh nodes

            % store variables
            nodes = self.mesh.nodes;
            if isempty(self.time)
                times = 0;
            else
                times = self.time.tGrid;
            end
            vals = zeros(length(nodes),length(times));

            % wrap function and evaluate at nodes
            func = Domain1d.wrapFunction(func);
            vals = func(nodes,times);
        end

        function vals = function2centroids(self,func)
        % FUNCTION2CENTROIDS Evaluate a function at the centroids of each element.
        % INPUTS:
        %   func - function to be evaluated (can be function handle, symfun,
        %           sym, or constant i.e. double) 
        % OUTPUTS:
        %   vals - function values at mesh nodes
        
            % store variables
            mids = self.mesh.centroids;
            if isempty(self.time)
                times = 0;
            else
                times = self.time.tGrid;
            end
            vals = zeros(length(mids),length(times));

            % wrap function and evaluate at centroids
            func = Domain1d.wrapFunction(func);
            vals = func(mids);

        end

        function vals = function2midpoints(self,func)
        % FUNCTION2MIDPOINTS Evaluate a function at the midpoints of each element.
        % INPUTS:
        %   func - function to be evaluated (can be function handle, symfun,
        %           sym, or constant i.e. double) 
        % OUTPUTS:
        %   vals - function values at mesh nodes
        %
        % NOTES:
        %   (1) In 1D, the centroids are the same as the midpoints. In 2d, i.e.
        %   in the Domain2d object, "centroids" are the centeroids of triangular
        %   elements while "midpoints" are the midpoints of the edges of the
        %   triangular elements. If we maintain this convention in 1D, then the
        %   midpoint of a cell should be referred to as its centroid and, as
        %   there are edges in 1D, there should be no midpoints. 
        %   So, technically, the function2centroids method is the "correct" name
        %   for this method in 1D, and functon2midpoints shouldn't exist, as
        %   there are no "midpoints" in 1D. However, I suspect that a user of
        %   the 1D code will think to use function2midpoints when they mean
        %   function2centroids, so both methods are provided for convenience.

            vals = self.function2centroids(func);

        end


        function vals = nodes2centroids(self,nodalVals)
        % NODES2CENTROIDS Interpolate nodal values to element centroids.
        % INPUTS:
        %   nodalVals - values at mesh nodes
        % OUTPUTS:
        %   vals      - interpolated values at element centroids

            % vectorize the above loop
            vals = 0.5 * (nodalVals(1:end-1) + nodalVals(2:end));

        end

        function vals = nodes2midpoints(self,nodalVals)
        % NODES2MIDPOINTS Interpolate nodal values to element midpoints.
        % INPUTS:
        %   nodalVals - values at mesh nodes
        % OUTPUTS:
        %   vals      - interpolated values at element midpoints
        %
        % NOTES:
        %   (1) See Notes (1) in function2midpoints method.

            vals = self.nodes2centroids(nodalVals);

        end


        % QUADRATURE
        function vals = twoPointQuadrature(self,func)
        % TWOPOINTQUADRATURE implements 2-point Gauss-Legendre quadrature on the
        %   stored 1D mesh.
        %
        % Inputs:
        %   func - function to be integrated (can be function handle, symfun,
        %          sym, or constant i.e. length 1 double, or vector of doubles)
        % Outputs:
        %   vals - integral values over each element

            % process input function
            func = self.twoPointQuadrature_evaluateFunction(func);
            fq = func.values;

            % get quadrature points and weights
            [xq, wq] = self.twoPointQuadrature_pointsAndWeights();

            % reshape fq into nElems x 2: each row = [f(q1), f(q2)] for that element
            fq_mat = reshape(fq, 2, []).';   % nElems x 2

            % elementwise integrals: sum_j wq(e,j) * f(e,j)
            elemVals = sum(fq_mat .* wq, 2); % nElems x 1

            % global integral = sum over elements
            vals = sum(elemVals);

        end

        function [xq, wq] = twoPointQuadrature_pointsAndWeights(self)
        % TWOPOINTQUADRATURE_POINTSANDWEIGHTS computes the quadrature points and
        %   weights for the stored mesh
        
            % store variables
            nodes  = self.mesh.nodes(:);
            elems  = self.mesh.elements;
            nElems = size(elems, 1);
            mids   = self.mesh.centroids;
            areas  = self.mesh.areas;

            % reference Gauss points and weights on [-1,1]
            xi_ref = [-1/sqrt(3); 1/sqrt(3)];
            w_ref  = [1; 1];

            % element endpoints
            x1 = nodes(elems(:,1));
            x2 = nodes(elems(:,2));

            % quad points and weights
            xq = [mids - 0.5*areas/sqrt(3), mids + 0.5*areas/sqrt(3)];
            wq = [0.5*areas, 0.5*areas];                  

        end

        function f_q = twoPointQuadrature_nodes2quadpoints(self,f_nodal)
        % TWOPOINTQUADRATURE_NODES2QUADPOINTS Interpolates nodal values to
        %   the 2-point quadrature points on each element.
        %
        % Inputs:
        %   f_nodal - NodalVec object containing values at mesh nodes,
        %               representing P1 function 
        % Outputs:
        %   f_q     - QuadVec object containing function values at quadrature points

            % process input
            validateattributes(f_nodal, {'NodalVec'}, {'scalar'});
            f_nodal = f_nodal.values;

            % store variables
            coords = self.mesh.nodes(:);    % nNodes x 1
            elems  = self.mesh.elements;    % nElems x 2
            areas  = self.mesh.areas;       % nElems x 1

            % element endpoints
            x1 = coords(elems(:,1));        % nElems x 1
            x2 = coords(elems(:,2));        % nElems x 1

            % nodal values of function
            f1 = f_nodal(elems(:,1));       % nElems x 1
            f2 = f_nodal(elems(:,2));       % nElems x 1

            % Evaluate shape functions at quadrature points (nElems x nQ)
            [xq, ~] = self.twoPointQuadrature_pointsAndWeights();
            phi1 = (x2 - xq) ./ areas;      % φ1(x_q)
            phi2 = (xq - x1) ./ areas;      % φ2(x_q)

            % Combine to get function values
            f_q = phi1 .* f1 + phi2 .* f2;  % nElems x nQ

            % reshape to column vector
            f_q = reshape(f_q, [], 1);    % (nElems*nQ) x 1
            f_q = QuadVec(f_q);

        end


        function f_q = twoPointQuadrature_evaluateFunction(self,func)
        % TWOPOINTQUADRATURE_EVALUATEFUNCTION Evaluates a function at the
        %   2-point quadrature points on each element.
        % INPUTS:
        %   func - function to be evaluated (can be function handle, symfun,
        %           sym, or constant i.e. double) 
        % OUTPUTS:
        %   f_q  - function values at quadrature points


            % if func is a QuadVec, assign directly
            if isa(func,'QuadVec')
                f_q = func;

            % if func is a NodalVec, interpolate to quad points
            elseif isa(func,'NodalVec')
                f_q = self.twoPointQuadrature_nodes2quadpoints(func);

            % if func is a function handle, evaluate at quad points
            elseif isa(func,'function_handle')

                % get quadrature points
                [xq, ~] = self.twoPointQuadrature_pointsAndWeights();

                % evaluate function at quad points
                func = Domain1d.wrapFunction(func);
                f_q = func(xq(:),0);  % time variable is dummy here
                f_q = QuadVec(f_q);

            % if func is symfun or sym, wrap and evaluate at quad points
            elseif isa(func,'symfun') || isa(func,'sym')

                % get quadrature points
                [xq, ~] = self.twoPointQuadrature_pointsAndWeights();

                % evaluate function at quad points
                func = Domain1d.wrapFunction(func);
                f_q = func(xq(:),0);
                f_q = QuadVec(f_q);

            % if func is a double of length 1, wrap and evaluate at quad points
            elseif isa(func,'double') && length(func) == 1

                % get quadrature points
                [xq, ~] = self.twoPointQuadrature_pointsAndWeights();

                % evaluate function at quad points
                func = Domain1d.wrapFunction(func);
                f_q = func(xq(:),0);
                f_q = QuadVec(f_q);

            % if func is double vector, see if its length matches nQuadPoints
            elseif isa(func,'double') && length(func) > 1
                
                % get quadrature points
                [xq, ~] = self.twoPointQuadrature_pointsAndWeights();

                % evaluate function at quad points
                ptCount = size(func,1);
                matchesQuads = (ptCount == length(xq));
                matchesNodes = (ptCount == self.mesh.nNodes);

                errMsg = ['Unable to confidently determine whether input is ' ...
                        'a double vector of values at nodes or at quadrature points. ' ...
                        'Please provide either a QuadVec or NodalVec object ' ...
                        'to avoid ambiguity.'];

                % if matches both or neither, throw error
                if matchesQuads && matchesNodes
                    error(errMsg);

                % if matches quadrature points, create QuadVec
                elseif matchesQuads
                    f_q = QuadVec(func);

                % if matches nodes, create NodalVec and interpolate to quad points
                elseif matchesNodes
                    nodalVec = NodalVec(func);
                    f_q = self.twoPointQuadrature_nodes2quadpoints(nodalVec);
                else
                    error(errMsg);
                end

            else
                error('Unsupported function type.');
            end

        end

        function vals = IP_twoPointQuadrature(self,func1,func2)
        % IP_TWOPOINTQUADRATURE Computes the inner product of two functions
        %   using 2-point Gauss-Legendre quadrature on the stored 1D mesh.
        % INPUTS:
        %   func1 - first function (can be function handle, symfun,
        %            sym, or constant i.e. double) 
        %   func2 - second function (can be function handle, symfun,
        %            sym, or constant i.e. double) 
        % OUTPUTS:
        %   vals  - inner product values

            % process input functions
            f1q = self.twoPointQuadrature_evaluateFunction(func1);
            f2q = self.twoPointQuadrature_evaluateFunction(func2);
            f1q = f1q.values;
            f2q = f2q.values;

            % take product at quadrature points
            fq = f1q .* f2q;
            fq = QuadVec(fq);

            % call quadrature method
            vals = self.twoPointQuadrature(fq);

        end

        function vals = L2norm_twoPointQuadrature(self,func)
        % L2NORM_TWOPOINTQUADRATURE Computes the L2 norm of a function
        %   using 2-point Gauss-Legendre quadrature on the stored 1D mesh.
        % INPUTS:
        %   func - function (can be function handle, symfun,
        %            sym, or constant i.e. double) 
        % OUTPUTS:
        %   vals - L2 norm value    

            % process input function
            fq = self.twoPointQuadrature_evaluateFunction(func);
            
            % call the inner product
            fq = IP_twoPointQuadrature(self, fq, fq);
            
            % take square root to get L2 norm
            vals = sqrt(fq);

        end

        function vals = L2err_twoPointQuadrature(self,func1,func2)
        % L2ERR_TWOPOINTQUADRATURE Computes the L2 error between two functions
        %   using 2-point Gauss-Legendre quadrature on the stored 1D mesh.
        % INPUTS:
        %   func1 - first function (can be function handle, symfun,
        %            sym, or constant i.e. double) 
        %   func2 - second function (can be function handle, symfun,
        %            sym, or constant i.e. double) 
        % OUTPUTS:
        %   vals  - L2 error value

            % process input functions
            f1q = self.twoPointQuadrature_evaluateFunction(func1);
            f2q = self.twoPointQuadrature_evaluateFunction(func2);

            % take difference at quadrature points
            fq = f1q.values - f2q.values;
            fq = QuadVec(fq);
            
            % call the L2 norm
            vals = L2norm_twoPointQuadrature(self, fq);

        end
            
    end

    methods (Static)        
        function func = wrapFunction(func)
        % WRAPFUNCTION Ensure that the input function is a function handle.
        % USAGE:
        %   func = wrapFunction(self,func)
        % INPUTS:
        %   func - function to be wrapped (can be function handle, symfun, sym, or constant)
        % OUTPUTS:
        %   func - wrapped function handle  
        %
        % NOTE:
        %   (1) Wrapping the function is necessary in the case where the input is
        %   constant either in time or in space. The expected output of any
        %   function evaluation is a matrix of size nSpaceEvals x nTimeEvals.
        %   But if a function is constant in space or time, MATLAB will return a
        %   single value an array of constants. Wrapping the function ensures that
        %   the output is always of the correct size. 
        %
        %   (2) A function is `wrapped' by adding an array of zeros of the same
        %   size as the expected output. However, the exact dimensions of this
        %   output array do not need to be known at the time of wrapping.
        %   Rather, by adding: 
        %       
        %       func = @(x,t) func(x,t) + zeros(x,t);
        %   
        %   a zero array of the correct size will be created at the time of
        %   function evaluation. This is sufficient to ensure that the output is
        %   always of the correct size.
        %
        %   (3) The expected input of any function evaluation is a
        %   function_handle, so wrapFunction also converts symfun, sym, or
        %   constant (i.e. double) inputs to a function handle. In the case of a
        %   constant function, evauating a function_handle is slightly more
        %   expensive than just creating an array of the appropriate size.
        %   However, this extra computational cost is worth it for the sake of
        %   reducing code complexity. 
        
            % if func is a constant (i.e. double) convert to function_handle
            if isa(func,'double') && length(func) == 1
                func = @(x,t) func;
            elseif isa(func,'double') && length(func) > 1
                error('If func is a double, it must be a scalar value.');
            end

            % if func is symfun convert to function_handle
            if isa(func,'symfun')
                syms x t
                func = matlabFunction(symfun(func,[x t]));
            end

            % if func is sym convert to function_handle
            if isa(func,'sym')
                syms x t
                func = matlabFunction(symfun(func,[x t]));
            end

            % if func is function_handle but only of x, convert to function_handle of (x,t)
            if isa(func,'function_handle') 
                func = Domain1d.checkVariables(func);
            end

            % wrap function handle
            func = @(x,t) func(x,t) + zeros(length(x),length(t));

        end

        function func = checkVariables(func)

            % store function type data
            isSpaceVarying = num2str(Coefficients.hasVariable(func,'x'));
            isTimeVarying  = num2str(Coefficients.hasVariable(func,'t'));
            code = strcat(isSpaceVarying, isTimeVarying);
            code = str2num(code);

            % update function signature based on variable dependence
            switch code

                % no variables
                case 00
                    func = @(x,t) func();
                
                % time variable only
                case 01
                    func = @(x,t) func(t);

                % space variable only
                case 10
                    func = @(x,t) func(x);

                % space and time variables
                case 11
                    func = @(x,t) func(x,t);
            end
        end


    end

end