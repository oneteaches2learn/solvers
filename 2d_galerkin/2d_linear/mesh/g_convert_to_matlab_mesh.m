clear all
% STEP 1: run the mesh as a .m script
%   the mesh will be automatically stored as msh
%   msh is a structure with mesh data, but not a mesh object
domain_dist100;


% STEP 2: store the elements; 
%   trim off unnecessary fourth column (bc elements are 3d)
%   transpose because matlab expects a 3 x nElem matrix
%   reorder the nodes to match matlab's ordering
elements = msh.TRIANGLES(:,1:3);
elements = elements';
elements = elements([1 3 2],:);


% STEP 3: store the nodes;
%   trim off the third column (bc nodes are 2d)
%   transpose because matlab expects a 2 x nNode matrix
nodes = msh.POS(:,1:2);
nodes = nodes';


% STEP 4: create the mesh object
%   to use matlab's pde toolbox, define a model
%   the mesh and geometry will be stored in the model
model = createpde();
geo = geometryFromMesh(model,nodes,elements);


% STEP 5: create Tyler's mesh object
mesh = Mesh2d();
mesh.base = 2;
mesh.p = 4;
mesh.h = 1 / mesh.base^mesh.p;
mesh.Mesh = model.Mesh;
mesh.effectiveNodes = [1:1:mesh.nNodes];
mesh.effectiveElems = [1:1:mesh.nElems];


% STEP 6: create Tyler's boundary object
bdry = Boundary2d();

% get endpoints of edge 20 from model.geometry
edges = zeros(model.Geometry.NumEdges,4);
for i = 1:model.Geometry.NumEdges
    nodes_i = findNodes(model.Mesh,"region","Edge",i);
    vert_1 = model.Mesh.Nodes(:,nodes_i(1))';
    vert_2 = model.Mesh.Nodes(:,nodes_i(end))';
    edges(i,:) = [vert_1, vert_2];
end
[edges,edgeDict] = orderEdges(edges);

% create edge list for boundary object
%{
% create one edge object for each edge
edgesObj = [];
for i = 1:size(edges,1)
    vert1 = edges(i,1:2);
    vert2 = edges(i,3:4);
    edge_i = BoundaryEdge2d(vert1,vert2);
    edge_i.segIDs = i;
    edge_i.nodes = findNodes(model.Mesh,"region","Edge",edgeDict(i));
    edgesObj = [edgesObj, edge_i];
end
%}

% hack: create four total edges
edgesObj = [];

% first three edges are the first three edges
for i = 1:3
    vert1 = edges(i,1:2);
    vert2 = edges(i,3:4);
    edge_i = BoundaryEdge2d(vert1,vert2);
    edge_i.segIDs = i;
    edge_i.nodes = findNodes(model.Mesh,"region","Edge",edgeDict(i));
    edgesObj = [edgesObj, edge_i];
end

% last edge contains all data from all remaining edges
vert1 = edges(4,1:2);
vert2 = edges(end,3:4);
edge_i = BoundaryEdge2d(vert1,vert2);
edge_i.segIDs = [4:1:size(edges,1)];
edge_i.nodes = [];
for i = 4:size(edges,1)
    edge_i.nodes = [edge_i.nodes, findNodes(model.Mesh,"region","Edge",edgeDict(i))];
end
edge_i.nodes = unique(edge_i.nodes);
edgesObj = [edgesObj, edge_i];

% add edges to boundary object
bdry = Boundary2d();
bdry.edges = edgesObj;


% STEP 7: create Tyler's domain object
dom = Domain2d();
dom.mesh = mesh;
dom.boundary = bdry;


function [ordered_edges,edgeDict] = orderEdges(edges)
% Reorder edges to form a counterclockwise cycle around the boundary
%
% INPUTS:
%   edges: n x 4 matrix, where each row is an edge defined by its two
%          endpoints
%
% OUTPUTS:
%   ordered_edges: m x 4 matrix, where each row is an edge defined by its
%                  two endpoints, and the edges form a counterclockwise cycle
%                  around the boundary
%   edgeDict: m x 1 vector, where each element is the index of the edge in
%             the original 'edges' matrix
    
    % STEP 1: find the starting vertex
    % extract vertex list from edge list
    vertices = unique([edges(:,1:2); edges(:,3:4)],'rows');
    
    % find the starting vertex (minimum y, then minimum x)
    min_y = min(vertices(:,2));
    min_y_points = vertices(vertices(:,2)==min_y,:);
    [min_x, idx_min_x] = min(min_y_points(:,1));
    startVertex = min_y_points(idx_min_x,:);
    

    % STEP 2: find the starting edge
    % find edges containing the starting vertex
    startEdges = find( (edges(:,1) == startVertex(1) & edges(:,2) == startVertex(2)) | ...
                       (edges(:,3) == startVertex(1) & edges(:,4) == startVertex(2)) );
                            
    % choose the edge representing counterclockwise traversal
    angles = zeros(length(startEdges),1);
    for i = 1:length(startEdges)

        % put startVertex at the begining of the edge
        idx = startEdges(i);
        edge = edges(idx,:);
        if edge(1:2) == startVertex
            otherVertex = edge(3:4);
        else
            otherVertex = edge(1:2);
            edges(idx,:) = [startVertex, otherVertex]; % Reverse edge
        end

        % calculate the angles these edges make with the x-axis
        delta = otherVertex - startVertex;
        angles(i) = atan2(delta(2), delta(1));
    end

    % choose the edge with the smallest angle to be the starting edge
    [~, minAngleIdx] = min(angles);
    startEdgeIdx = startEdges(minAngleIdx);
    startEdge = edges(startEdgeIdx,:);

    % update the ordered list and dictionary; mark the starting edge as used
    ordered_edges = startEdge;
    used_edges = false(size(edges,1),1);
    used_edges(startEdgeIdx) = true;
    edgeDict = [startEdgeIdx];
    
    
    % STEP 3: traverse the remaining edges
    currentVertex = startEdge(3:4);
    while sum(used_edges) < size(edges,1)

        % find the remaining edge that is connected to the current vertex
        connectedEdgeIdx = find(~used_edges & (...
            (edges(:,1) == currentVertex(1) & edges(:,2) == currentVertex(2)) | ...
            (edges(:,3) == currentVertex(1) & edges(:,4) == currentVertex(2))));
        if isempty(connectedEdgeIdx)
            break
        end

        % put vertices of connected edge in the correct order
        connectedEdge = edges(connectedEdgeIdx,:);
        if connectedEdge(1) == currentVertex(1) && connectedEdge(2) == currentVertex(2)
            nextVertex = connectedEdge(3:4);
        else
            nextVertex = connectedEdge(1:2);
            connectedEdge = [currentVertex, nextVertex]; % Reverse edge
        end

        % store the connected edge and mark it as used
        ordered_edges = [ordered_edges; connectedEdge];
        used_edges(connectedEdgeIdx) = true;
        edgeDict = [edgeDict; connectedEdgeIdx];

        % update the current vertex
        currentVertex = nextVertex;

    end

    % The 'ordered_edges' variable now contains the edges in counterclockwise order
end
    