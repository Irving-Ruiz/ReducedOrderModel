function u_vol = volume_average(nodes, elements, u)
% VOLUME_AVERAGE: Computes the volume-averaged value of a scalar field
% over a tetrahedral mesh.
%
% Input arguments:
% nodes: Node coordinates (3 x N)
% elements: Element connectivity (4 x Ne)
% u: Scalar field evaluated at the nodes (1 x N or compatible indexing)
%
% Output arguments:
% u_vol: Volume-averaged value of the scalar field

% Number of elements (not explicitly used, but kept for clarity)
E = size(elements,2);

% Coordinates of the four nodes of each tetrahedron
X1 = nodes(:,elements(1,:));
X2 = nodes(:,elements(2,:));
X3 = nodes(:,elements(3,:));
X4 = nodes(:,elements(4,:));

% Compute tetrahedral volumes using scalar triple product
V = abs(dot(X2 - X1, cross(X3 - X1, X4 - X1))) / 6;

% Total volume of the domain
Vtotal = sum(V);

% Compute element-averaged values (linear interpolation from nodal values)
Ue = (u(:,elements(1,:)) + u(:,elements(2,:)) + ...
      u(:,elements(3,:)) + u(:,elements(4,:))) / 4;

% Compute volume-weighted average
u_vol = (Ue * V') / Vtotal;

end