function dudt = mass_transfer_equilibrium_model(t, u, nNodes, ...
    interior_nodes, a, K, itersol, nodes, elements, boundary_nodes, ...
    stencil)
% MASS_TRANSFER_EQUILIBRIUM_MODEL: Computes time derivative of concentration
% field considering diffusion and equilibrium with external phase.
%
% Input arguments:
% t: Time
% u: Solution vector at interior nodes
% nNodes: Total number of nodes
% interior_nodes: Indices of interior nodes
% a: Solution-to-solid volume ratio
% K: Partition coefficient
% itersol: Maximum iterations for equilibrium update
% nodes, elements: Mesh definition
% boundary_nodes: Indices of boundary nodes
% stencil: Precomputed stencil structure
%
% Output:
% dudt: Time derivative at interior nodes

% Number of interior nodes
nInterior = length(interior_nodes);

% Initialize derivative
dudt = zeros(nInterior,1);

% Reconstruct full solution (interior + boundary)
ufull = zeros(nNodes,1);
ufull(interior_nodes) = u(:);
ufull(boundary_nodes) = 0;

% Equilibrium coupling with external phase
if isinf(a)
    % Infinite external volume → no iteration required
    u_vol = volume_average(nodes,elements,ufull');
else
    % Finite volume → iterative update of fluid concentration
    if t == 0
        CL = 0; % Initial fluid concentration
    else
        CL = 0;
        tol = 1e-8; % Convergence tolerance

        for i = 1:itersol
            % Compute volume-averaged concentration
            u_vol = volume_average(nodes,elements,ufull');

            % Update fluid concentration
            CL_new = K*(1 - u_vol)/a;

            % Convergence check
            if abs(CL_new - CL) < tol
                break
            end

            % Update boundary condition
            CL = CL_new;
            ufull(boundary_nodes) = CL/K - 1/a;

            % Debug message if max iterations reached
            if i == itersol
                "llegue al lím"
            end
        end
    end
end

% Compute spatial derivatives using GFD stencil
derivs = compute_derivatives_equilibrium(ufull',stencil);

% Laplacian (sum of second derivatives)
dudt(:) = derivs(4,interior_nodes) + ...
          derivs(5,interior_nodes) + ...
          derivs(6,interior_nodes);

end