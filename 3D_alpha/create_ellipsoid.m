function [nodes, elements, nNodes, nElements, status] = ...
    create_ellipsoid(a, b, c, hmax)
% CREATE_ELLIPSOID: Generates a tetrahedral mesh of an ellipsoid by scaling
% a unit sphere.
%
% Input arguments:
% a, b, c: Semi-axes of the ellipsoid
% hmax: Maximum element size (mesh resolution parameter)
%
% Output arguments:
% nodes: Node coordinates (3 x N)
% elements: Element connectivity
% nNodes: Number of nodes
% nElements: Number of elements
% status: Structure with mesh information (from node cleanup)

% Create PDE model
model = createpde;

% Generate unit sphere geometry
gm = multisphere(1);  % radius = 1

% Scale sphere to ellipsoid with semi-axes (a, b, c)
scale(gm,[a b c]);

% Assign geometry to model
model.Geometry = gm;

% Generate tetrahedral mesh with prescribed maximum element size
generateMesh(model,Hmax=hmax);

% Extract mesh data
nodes = model.Mesh.Nodes;
elements = model.Mesh.Elements;

% Remove unused nodes and update connectivity
[nodes, elements, status] = remove_unused_nodes(nodes, elements);

% Extract final mesh sizes
nNodes = status.nNodes;
nElements = status.nElements;

end