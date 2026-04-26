function stencil = precompute_stencils(nodes,first_level,second_level)
% PRECOMPUTE_STENCILS: Builds local GFD stencils for each node
%
% Input arguments:
% nodes: Node coordinates (3 x N)
% first_level: First-level neighbors for each node
% second_level: Second-level neighbors for each node
%
% Output arguments:
% stencil: Structure containing, for each node:
%   - neigh: Neighbor indices
%   - w: Weights
%   - A_mat: Local system matrix
%   - L_mat, U_mat: LU factors of A_mat
%   - cond_A: Condition number of A_mat

% Threshold for acceptable conditioning of local system
cond_A_tol = 1E12;

% Number of nodes
N = size(nodes,2);

% Initialize stencil structure
stencil = struct( ...
    'neigh',cell(N,1), ...
    'w',cell(N,1), ...
    'A_mat',cell(N,1), ...
    'L_mat',cell(N,1), ...
    'U_mat',cell(N,1), ...
    'cond_A',zeros(N,1) );

% Loop over all nodes
for k = 1:N
    
    % Initial neighbor set (first level)
    neigh = first_level{k};

    % Ensure sufficient neighbors for second-order reconstruction
    if length(neigh) < 9
        neigh = unique([first_level{k} second_level{k}]);
    end

    % Compute geometric differences and weights
    [dx, dy, dz, ~, w] = compute_weighs(nodes, neigh, k);

    % Assemble local system matrix
    [A, cond_A] = compute_matrix(dx, dy, dz, w);

    % If ill-conditioned, expand stencil and recompute
    if cond_A > cond_A_tol
        neigh = unique([first_level{k} second_level{k}]);
        [dx, dy, dz, ~, w] = compute_weighs(nodes, neigh, k);
        [A, cond_A] = compute_matrix(dx, dy, dz, w);
    end

    % LU factorization for efficient reuse
    [L_mat,U_mat] = lu(A);

    % Store stencil data
    stencil(k).neigh = neigh;
    stencil(k).w     = w;
    stencil(k).A_mat = A;
    stencil(k).cond_A = cond_A;
    stencil(k).L_mat = L_mat;
    stencil(k).U_mat = U_mat;
end

end

%% Compute weighs
function [dx, dy, dz, r, w] = compute_weighs(nodes, neigh, k)
% COMPUTE_WEIGHS: Computes relative coordinates and weights for neighbors

% Relative coordinates with respect to node k
dx = nodes(1,neigh) - nodes(1,k);
dy = nodes(2,neigh) - nodes(2,k);
dz = nodes(3,neigh) - nodes(3,k);

% Euclidean distance
r = sqrt(dx.^2 + dy.^2 + dz.^2);

% Weight function (inverse distance squared)
% Alternative options (commented):
% w = 1 ./ (r.^3); % can be unstable
% w = exp(-(r./max(r)).^2);
w = 1 ./ (r.^2);

end

%% Matrix components
function [A, cond_A] = compute_matrix(dx, dy, dz, w)
% COMPUTE_MATRIX: Builds local system matrix for GFD approximation

% First-order terms
A1 = (w.*dx)';
A2 = (w.*dy)';
A3 = (w.*dz)';

% Second-order terms (diagonal)
A4 = (0.5*w.*dx.^2)';
A5 = (0.5*w.*dy.^2)';
A6 = (0.5*w.*dz.^2)';

% Mixed second-order terms
A7 = (w.*dx.*dy)';
A8 = (w.*dx.*dz)';
A9 = (w.*dy.*dz)';

% Assemble matrix
A = [A1 A2 A3 A4 A5 A6 A7 A8 A9];

% Condition number (measure of numerical stability)
cond_A = cond(A);

end