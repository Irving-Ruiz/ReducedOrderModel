function [first_level, second_level] = ...
    build_node_neighbor_levels(elements, nNodes)
% BUILD_NODE_NEIGHBOR_LEVELS: Constructs first- and second-level neighbor
% lists for each node based on mesh connectivity.
%
% Input arguments:
% elements: Element connectivity (4 x Ne)
% nNodes: Total number of nodes
%
% Output arguments:
% first_level: Cell array where each entry contains first-level neighbors
% second_level: Cell array where each entry contains second-level neighbors

% Number of elements
Ne = size(elements,2);

% Initialize first-level neighbors
first_level = cell(nNodes,1);

% Loop over all elements
for e = 1:Ne
    
    % Nodes of current tetrahedron
    n = elements(:,e);
    
    % Each node is connected to the other nodes in the element
    for i = 1:4
        for j = 1:4
            if i ~= j
                first_level{n(i)} = [first_level{n(i)} n(j)];
            end
        end
    end
end

% Remove duplicate neighbors
for i = 1:nNodes
    first_level{i} = unique(first_level{i});
end

% Initialize second-level neighbors
second_level = cell(nNodes,1);

% Build second-level neighbors
for i = 1:nNodes

    % First-level neighbors of node i
    n1 = first_level{i};
    
    % Collect neighbors of neighbors
    n2 = [];

    for k = n1
        n2 = [n2 first_level{k}];
    end

    % Remove duplicates
    n2 = unique(n2);

    % Remove self node
    n2(n2 == i) = [];

    % Remove first-level neighbors
    n2 = setdiff(n2,n1);

    % Store second-level neighbors
    second_level{i} = n2;
end

end