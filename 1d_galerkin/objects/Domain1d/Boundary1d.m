classdef Boundary1d

    properties
        nEdges
        edges = BoundaryEdge1d.empty(1,0); % array of BoundaryEdge1d objects
        D_nodes
        N_nodes
        R_nodes
        P_nodes
        freeNodes
    end

    properties (Dependent)
        bTypes
        bConditions
        bConditions_du
    end

    properties (Hidden)
        nodes
    end

    methods
        % CONSTRUCTOR
        function self = Boundary1d(xLim, bcTypes, bcVals)

            if nargin > 0

                % preallocate edge array
                self.edges(2) = BoundaryEdge1d;

                % set edges
                self.edges(1) = BoundaryEdge1d(xLim(1));
                self.edges(2) = BoundaryEdge1d(xLim(2));

                % set outward normals
                self.edges(1).outwardNormal = -1;
                self.edges(2).outwardNormal = 1;

            end

            % if bcTypes provided, set them
            if nargin >= 2
                self.bTypes = bcTypes;
            end

            % if bcVals provided, set them
            if nargin == 3
                self.bConditions = bcVals;
            end

        end

        % SETTERS
        function self = set.bTypes(self,bcTypes)

            % set boundaryType on each edge
            for i = 1:self.nEdges
                self.edges(i).boundaryType = bcTypes(i);
            end
        end

        function self = set.bConditions(self,BCs)

            % validate BCs
            validateattributes(BCs, {'cell'}, {'vector'});   % must be a cell array (vector)
            assert(numel(BCs) == 2, 'Input must be a cell array of length 2.');

            % assign BCs
            for i = 1:length(BCs)
                self.edges(i).boundaryCondition = BCs{i};
            end

        end

        function self = set.bConditions_du(self,BCs_du)

            % validate BCs
            validateattributes(BCs_du, {'cell'}, {'vector'});   % must be a cell array (vector)
            assert(numel(BCs_du) == 2, 'Input must be a cell array of length 2.');

            % assign BCs
            for i = 1:length(BCs_du)
                self.edges(i).boundaryCondition_ddu = BCs_du{i};
            end

        end


        % GETTERS
        function val = get.nEdges(self)
            val = length(self.edges);
        end

        function val = get.bTypes(self)

            % Collect boundaryType from each edge
            if isempty(self.edges)
                val = '';
                return;
            end

            % self.edges(k).boundaryType is a char
            bt = arrayfun(@(e)e.boundaryType, self.edges, 'UniformOutput', false);
            val = [bt{:}];   % combine to e.g. 'DR'
        end

        function val = get.bConditions(self)
            
            % Collect boundaryCondition from each edge
            if isempty(self.edges)
                val = {};
                return;
            end

            % self.edges(k).boundaryCondition can be function handle, numeric, or cell
            bc = arrayfun(@(e)e.boundaryCondition, self.edges, 'UniformOutput', false);
            val = bc;   % combine to cell array
        end

        function val = get.bConditions_du(self)
            
            % Collect boundaryCondition from each edge
            if isempty(self.edges)
                val = {};
                return;
            end

            % self.edges(k).boundaryCondition can be function handle, numeric, or cell
            bc = arrayfun(@(e)e.boundaryCondition_ddu, self.edges, 'UniformOutput', false);
            val = bc;   % combine to cell array
        end


        % UTILITY METHODS
        function self = mesh2nodes(self,mesh)

            % assign nodes numbers from mesh to each edge
            self.edges(1).nodes = 1;
            self.edges(2).nodes = length(mesh.nodes);

            % clear out node lists
            self.D_nodes = [];
            self.N_nodes = [];
            self.R_nodes = [];
            self.P_nodes = [];

            % loop over edges and assign nodes based on boundary type
            for i = 1:self.nEdges
                edge = self.edges(i);
                switch edge.boundaryType
                    case 'D'
                        self.D_nodes = [self.D_nodes, edge.nodes];
                    case 'N'
                        self.N_nodes = [self.N_nodes, edge.nodes];
                    case 'R'
                        self.R_nodes = [self.R_nodes, edge.nodes];
                    case 'P'
                        self.P_nodes = [self.P_nodes, edge.nodes];
                end
            end

            % set free nodes
            boundNodes = [self.D_nodes, self.R_nodes];
            self.freeNodes = setdiff(mesh.effectiveNodes, boundNodes);

        end
            
    end
end


