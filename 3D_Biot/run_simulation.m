clc;
clear;
close all;

% Mass biot number
biot = 100;

% Parameter controlling mesh resolution (higher values increase resolution)
factorGeom = 11.5;

% Parameter controlling the simulation time scale (increase for lower Biot 
% values)
timeFactor = 1;

load('data.mat');

tsimLength = 400;

nGeom = length(dimensions);

figIndex = [1 19];
%figIndex = 1:nGeom;
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
    [model, nodes, elements, nNodes, nElements, nFaces] = ...
        create_ellipsoid(abc(1), abc(2), abc(3), hmax);

    % Specify boundary condition parameters
    if isinf(biot)
        applyBoundaryCondition(model,"dirichlet","Face",1:nFaces,"r",0);
    else
        applyBoundaryCondition(model,"neumann","Face",1:nFaces,"g",0, ...
            "q",biot);
    end

    % Specify model parameters
    specifyCoefficients(model,"m",0,"d",1,"c",1,"a",0,"f",0);

    % Specify initial condition
    setInitialConditions(model,1);

    % Final time
    tf = timeFactor*interp1([1;3],[4.75;1.5],sum(1./abc.^2));

    tsim = [linspace(0,tf,tsimLength)]';

    % Solve model
    results = solvepde(model,tsim);

    % Compute volume average solution
    up = volume_average(nodes,elements,results.NodalSolution');

    % Metrics
    uVolMin = min(up);

    % Store solution
    sol(j).res = [tsim, up];

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
save('results.mat',"res","sol","biot","factorGeom","timeFactor");