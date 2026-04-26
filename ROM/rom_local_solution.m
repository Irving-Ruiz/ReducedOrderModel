function y = rom_local_solution(nter, B, a, v, x0, xi, Fo, lambda_n)
% ROM_LOCAL_SOLUTION: Computes the local solution of the ROM
%
% Input arguments:
% nter: Number of eigenvalues (terms) to compute
% B: Biot number (dimensionless)
% a: Ratio of the volumes of solution and solid
% v: Bessel function order parameter
% x0: Initial guess for the first eigenvalue
% xi: Dimensionless position at which the solution is evaluated
% Fo: Fourier number (dimensionless time)
% lambda_n: Eigenvalues of the characteristic equation
%
% Output arguments:
% y: Solution evaluated at the specified xi and Fo

nt = length(Fo);           % Number of time points
y = zeros(nt,1);           % Preallocate solution vector

% Compute eigenvalues only if not provided
if isempty(lambda_n)
    lambda_n = rom_eigenvalues(nter, B, a, v, x0);
end

% Precompute coefficients
if isinf(B) && isinf(a)
    Cn = 2 ./ lambda_n ./ besselj(v+1, lambda_n);
    Xn = besselj(v, lambda_n * xi) ./ xi^v;

elseif isinf(B) && ~isinf(a)
    Cn = 2*(1+a) ./ ( ...
        a*besselj(v+1, lambda_n).*lambda_n + ...
        (2*v+2)^2*besselj(v+1, lambda_n)./lambda_n - ...
        (2*v+2)*besselj(v, lambda_n) );
    Xn = besselj(v, lambda_n * xi) ./ xi^v;

elseif ~isinf(B) && isinf(a)
    Cn = 2 .* besselj(v+1, lambda_n) ./ lambda_n ./ ...
        (besselj(v, lambda_n).^2 - ...
         besselj(v-1, lambda_n).*besselj(v+1, lambda_n));
    Xn = besselj(v, lambda_n * xi) ./ xi^v;

else
    error('Unsupported combination of parameters B and a.');
end

% Evaluate the series solution
for j = 1:nt
    if Fo(j) == 0
        y(j) = 1;
    else
        Tn = exp(-lambda_n.^2 * Fo(j));
        y(j) = sum(Cn .* Xn .* Tn);
    end
end

end