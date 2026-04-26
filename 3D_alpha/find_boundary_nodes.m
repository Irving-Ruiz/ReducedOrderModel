function boundary_nodes = find_boundary_nodes(elements)
% FIND_BOUNDARY_NODES: Identifies boundary nodes from tetrahedral mesh
% connectivity.
%
% Input arguments:
% elements: Element connectivity (4 x Ne)
%
% Output arguments:
% boundary_nodes: Indices of nodes located on the boundary

% Number of elements
Ne = size(elements,2);

% Each tetrahedron has 4 triangular faces
faces = zeros(4*Ne,3);
count = 0;

% Loop over elements to extract all faces
for e = 1:Ne
    
    % Nodes of current tetrahedron
    n = elements(:,e);
    
    % Define its 4 triangular faces
    tri = [n([1 2 3])'; 
           n([1 2 4])'; 
           n([1 3 4])'; 
           n([2 3 4])'];
    
    % Store sorted faces (for consistent comparison)
    for k = 1:4
        count = count + 1;
        faces(count,:) = sort(tri(k,:));
    end
end

% Sort faces to group identical ones
faces = sortrows(faces);

% Identify unique faces and their occurrences
[unique_faces,~,ic] = unique(faces,'rows');
counts = accumarray(ic,1);

% Boundary faces appear only once
boundary_faces = unique_faces(counts == 1,:);

% Extract boundary node indices
boundary_nodes = unique(boundary_faces(:));

end