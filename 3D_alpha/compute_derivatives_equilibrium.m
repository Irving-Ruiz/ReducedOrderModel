function derivs = compute_derivatives_equilibrium(u, stencil)
% COMPUTE_DERIVATIVES_EQUILIBRIUM: Computes spatial derivatives at each node
% using precomputed GFD stencils.
%
% Input arguments:
% u: Scalar field evaluated at nodes
% stencil: Structure containing local stencil data (neighbors, weights, 
% matrices)
%
% Output arguments:
% derivs: Matrix of derivatives (9 x N), where each column contains:
% [du/dx, du/dy, du/dz, d2u/dx2, d2u/dy2, d2u/dz2, d2u/dxdy, d2u/dxdz, 
% d2u/dydz]

% Number of nodes
N = length(stencil);

% Preallocate derivatives matrix
derivs = zeros(9,N);

% Loop over all nodes
for k = 1:N
    
    % Neighbor indices and weights
    neigh = stencil(k).neigh;
    w     = stencil(k).w;

    % Difference between neighbor values and central node
    du = u(neigh) - u(k);

    % Right-hand side of local system
    B = (w .* du)';

    % Solve local least-squares system to obtain derivatives
    % Alternative (using LU factors):
    % y = stencil(k).L_mat \ B;
    % X = stencil(k).U_mat \ y;
    
    X = stencil(k).A_mat \ B;

    % Store derivatives
    derivs(:,k) = X;
end

end