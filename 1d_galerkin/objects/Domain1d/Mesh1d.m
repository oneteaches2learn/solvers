classdef Mesh1d

	properties
		base           % mesh base, i.e. h = base^-p is max element size
		p              % mesh power, i.e. h = base^-p is max element size
		h              % mesh size, i.e. h = base^-p is max element size
		nNodes         % number of nodes
		nElems         % number of elements
		nodes          % node coordinates
		elements       % element connectivity matrix
		elementCoords  % element coordinates
        effectiveNodes % list of node indices used in effective domain
		centroids      % list of element midpoints
        midpoints      % list of element midpoints
		areas          % element areas
	end

    methods
        % CONSTRUCTOR
        function self = Mesh1d(x,base,p)

            % generate uniform mesh from base^-p
            if nargin == 3

                % store inputs
                self.p = p;
                self.base = base;
                
                % compute mesh size
                h = base^-p;
            
                % generate uniform mesh
                self.nodes = (x(1):h:x(2))';

            end
        end

        % GETTERS
        function val = get.nNodes(self)
            val = length(self.nodes);
        end

        function val = get.nElems(self)
            val = length(self.nodes) - 1;
        end

        function val = get.elements(self)
            val = [(1:self.nElems)', (2:self.nElems+1)'];
        end

        function val = get.elementCoords(self)
            val = [self.nodes(1:end-1), self.nodes(2:end)];
        end

        function val = get.effectiveNodes(self)
            val = 1:self.nNodes;
        end

        function val = get.centroids(self)
            val = 0.5 * (self.elementCoords(:, 1) + self.elementCoords(:, 2));
        end

        function val = get.midpoints(self)
            val = self.centroids;
        end

        function val = get.areas(self)
            val = self.nodes(2:end) - self.nodes(1:end-1);
        end

        function val = get.h(self)
            val = max(self.areas);
        end

    end
end
 