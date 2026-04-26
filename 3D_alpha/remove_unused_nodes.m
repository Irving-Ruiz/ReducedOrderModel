function [nodes, elements, status] = remove_unused_nodes(nodes, elements)
% REMOVE_UNUSED_NODES: Removes nodes not referenced by any element and
% updates mesh connectivity accordingly.
%
% Input arguments:
% nodes: Node coordinates (3 x N)
% elements: Element connectivity (at least 4 x Ne)
%
% Output arguments:
% nodes: Filtered node coordinates
% elements: Updated connectivity with reindexed nodes
% status: Structure with mesh information

% Original number of nodes
nNodes0 = size(nodes,2);

% Keep only tetrahedral connectivity (first 4 nodes per element)
elements = elements(1:4,:);

% Original number of elements (not explicitly used)
nElements0 = size(elements,2);

% Identify nodes actually used in the mesh
used_nodes = unique(elements(:));

% Create mapping from old node indices to new indices
map = zeros(1, size(nodes,2));
map(used_nodes) = 1:length(used_nodes);

% Remove unused nodes
nodes = nodes(:, used_nodes);

% Updated number of nodes
nNodes = size(nodes,2);

% Update element connectivity with new indexing
elements = map(elements);

% Updated number of elements
nElements = size(elements,2);

% Store status information
status.mesh_flag = (nNodes < nNodes0); % true if nodes were removed
status.nNodes0 = nNodes0;              % original number of nodes
status.nNodes = nNodes;                % updated number of nodes
status.nElements = nElements;          % number of elements

end