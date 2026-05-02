function [model, nodes, elements, nNodes, nElements, nFaces] = ...
    create_ellipsoid(a, b, c, hmax)
% CREATE_ELLIPSOID: Generates a tetrahedral mesh of an ellipsoid by scaling
% a unit sphere.
%
% Input arguments:
% a, b, c: Semi-axes of the ellipsoid
% hmax: Maximum element size (mesh resolution parameter)
%
% Output arguments:
% model: PDE model structure
% nodes: Node coordinates (3 x N)
% elements: Element connectivity
% nNodes: Number of nodes
% nElements: Number of elements
% nFaces: Number of faces

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
elements = elements(1:4,:);
nNodes = size(nodes,2);
nElements = size(elements,2);
nFaces = model.Geometry.NumFaces;

end