function [v] = gfd3d_ellipsoid(nodes, elements, a, K, ...
    itersol, tfinal, stencil)
% GFD3D_ELLIPSOID: Solves the transient mass transfer problem in
% an ellipsoidal domain using a GFD spatial discretization and ODE
% integration.
%
% Input arguments:
% nodes, elements: Mesh definition
% a: Solution-to-solid volume ratio
% K: Partition coefficient
% itersol: Maximum iterations for equilibrium update
% tfinal: Final simulation time
% stencil: Precomputed stencil structure for GFD derivatives
%
% Output arguments:
% v: Structure containing:
%   - t: Time vector
%   - u: Solution at interior nodes
%   - CL: Fluid-phase concentration
%   - up: Volume-averaged concentration
%   - uVolEq: Equilibrium concentration

% Total number of nodes
nNodes = size(nodes,2);

% Identify boundary nodes
boundary_nodes = find_boundary_nodes(elements);

% Logical mask for boundary nodes
isBoundary = false(1,nNodes);
isBoundary(boundary_nodes) = true;

% Interior nodes
interior_nodes = find(~isBoundary);
nInterior = length(interior_nodes);

% Initial condition (uniform concentration)
u0 = ones(1,nInterior);

% Time integration of semi-discrete system
[t,u] = ode113(@(t,u) mass_transfer_equilibrium_model(t, u, nNodes, ...
    interior_nodes, a, K, ...
    itersol, nodes, elements, boundary_nodes, stencil), ...
    [0 tfinal], u0);

% Post-processing: volume-averaged concentration and fluid phase
u_vol = zeros(length(t),1);
CL = zeros(length(t),1);

% Equilibrium volume-averaged concentration
uVolEq = 1/(1+a);

for j = 1:length(t)
    
    % Reconstruct full solution
    ufull = zeros(nNodes,1);
    ufull(interior_nodes) = u(j,:);
    ufull(boundary_nodes) = 0;

    if isinf(a)
        % Infinite external volume
        u_vol(j) = volume_average(nodes,elements,ufull');
    else
        % Finite volume: iterative equilibrium coupling
        tol = 1e-8;

        for i = 1:itersol
            
            % Compute volume-averaged concentration
            u_vol(j) = volume_average(nodes,elements,ufull');
            
            % Update fluid concentration
            CL_new = K*(1 - u_vol(j))/a;

            % Convergence check
            if abs(CL_new - CL(j)) < tol
                break
            end

            % Update boundary condition
            CL(j) = CL_new;
            ufull(boundary_nodes) = CL(j)/K - 1/a;
        end
    end
end

% Store results
v.t = t;         % Time
v.u = u;         % Interior solution
v.CL = CL;       % Fluid-phase concentration
v.up = u_vol;    % Volume-averaged concentration
v.uVolEq = uVolEq; % Equilibrium value

end