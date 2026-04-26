function [lambda_n, dlambda_n] = rom_eigenvalues(nter, B, a, v, x0)
% ROM_EIGENVALUES: Computes eigenvalues of the ROM
%
% Input arguments:
% nter: Number of eigenvalues (terms) to compute
% B: Biot number (dimensionless)
% a: Ratio of the volumes of solution and solid
% v: Bessel function order parameter
% x0: Initial guess for the first eigenvalue
%
% Output arguments:
% lambda_n: Vector of computed eigenvalues (roots of the characteristic equation)
% dlambda_n: Differences between consecutive eigenvalues (used for validation;
% spacing is expected to be approximately pi)

lambda_n = zeros(nter,1);  % Preallocate eigenvalues

for n = 1:nter
    if n == 1
        x_init = x0;
    else
        % Search interval for subsequent eigenvalues
        x_init = [lambda_n(n-1) + 2.5, lambda_n(n-1) + 3.5];
    end
    
    % Solve the characteristic equation
    lambda_n(n) = fzero(@(x) flambdan(x, B, a, v), x_init);
end

% Compute spacing between successive eigenvalues
dlambda_n = diff(lambda_n);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function f = flambdan(x, B, a, v)
% FLAMBDAN: Characteristic equation

if isinf(B) && isinf(a)
    f = besselj(v, x);

elseif isinf(B) && ~isinf(a)
    f = (a*x^2 + (2*v+2)*v) * besselj(v, x) - ...
        (2*v+2)*x*(besselj(v-1, x) - besselj(v+1, x))/2;

elseif ~isinf(B) && isinf(a)
    f = x*(besselj(v-1, x) - besselj(v+1, x))/2 + ...
        (B - v)*besselj(v, x);

else
    error('Unsupported combination of parameters B and a.');
end

end