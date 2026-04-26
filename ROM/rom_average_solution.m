function Y = rom_average_solution(nter, B, a, v, x0, Fo, lambda_n)
% ROM_AVERAGE_SOLUTION: Computes the volume-averaged solution of the ROM
%
% Input arguments:
% nter: Number of eigenvalues (terms) to compute
% B: Biot number (dimensionless)
% a: Ratio of the volumes of solution and solid
% v: Bessel function order parameter
% x0: Initial guess for the first eigenvalue
% Fo: Fourier number (dimensionless time)
% lambda_n: Eigenvalues of the characteristic equation
%
% Output arguments:
% Y: Volume-averaged solution evaluated at Fo

nt = length(Fo);           % Number of time points
Y = zeros(nt,1);           % Preallocate solution vector

% Compute eigenvalues only if not provided
if isempty(lambda_n)
    lambda_n = rom_eigenvalues(nter, B, a, v, x0);
end

% Precompute coefficients
if isinf(B) && isinf(a)
    An = 2 ./ lambda_n.^2;

elseif isinf(B) && ~isinf(a)
    nume = 2*(1+a);
    deno = a*besselj(v+1, lambda_n).*lambda_n + ...
           (2*v+2)^2*besselj(v+1, lambda_n)./lambda_n - ...
           (2*v+2)*besselj(v, lambda_n);
    An = nume ./ deno;
    An = An .* besselj(v+1, lambda_n) ./ lambda_n;

elseif ~isinf(B) && isinf(a)
    An = 2 .* besselj(v+1, lambda_n).^2 ./ lambda_n.^2 ./ ...
        (besselj(v, lambda_n).^2 - ...
         besselj(v-1, lambda_n).*besselj(v+1, lambda_n));

else
    error('Unsupported combination of parameters B and a.');
end

% Evaluate the series solution
for j = 1:nt
    if Fo(j) < 1e-12
        Y(j) = 1;
    else
        Tn = exp(-lambda_n.^2 * Fo(j));
        Y(j) = (2*v+2) * sum(An .* Tn);
    end
end

end