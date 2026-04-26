clc;
clear;
close all;

% Solution-to-solid volume ratio
alpha = inf;

% Partition coefficient
K = 1;

% Maximum number of iterations for updating the fluid-phase concentration
itersol = 20;

% Parameter controlling mesh resolution (higher values increase resolution)
factorGeom = 10;

% Parameter controlling the simulation time scale (reduce for lower alpha 
% values)
timeFactor = 1;

load('data.mat');

nGeom = size(dimensions,1);

figIndex = 1:nGeom;
nLoop = length(figIndex);

tic;

for j = 1:nLoop
    
    disp(j)
    t1 = tic;

    iFig = figIndex(j);

    % Normalize geometry
    abc = sort(dimensions(iFig,:) / min(dimensions(iFig,:)));

    % Mesh size
    hmax = geomean(abc) / factorGeom;

    % Mesh generation
    [nodes, elements, nNodes, nElements, status] = ...
        create_ellipsoid(abc(1), abc(2), abc(3), hmax);

    % Neighbor structure
    [first_level, second_level] = ...
        build_node_neighbor_levels(elements, nNodes);

    % Stencil
    stencil = precompute_stencils(nodes, first_level, second_level);

    % Final time
    tf = timeFactor * interp1([1;3], [4.75;1], sum(1./abc.^2), ...
        'linear', 'extrap');

    % Solve model
    results = gfd3d_ellipsoid(nodes, elements, alpha, K, ...
        itersol, tf, stencil);

    % Store solution
    sol(j).res = [results.t, results.up];

    % Metrics
    uVolMin = min(results.up);

    % Summary
    res(j,:) = [abc, hmax, tf, nNodes, nElements, uVolMin];

    time = toc(t1);
    disp([j time])

end

toc

% Plot results
for j = 1:nLoop
    figure;
    semilogy(sol(j).res(:,1), sol(j).res(:,2), 'LineWidth', 1.5);
    xlabel('Fo');
    ylabel('\Psi_{avg}');
    grid on;
end

% Save results
save('results.mat',"res","sol","alpha","factorGeom","timeFactor");