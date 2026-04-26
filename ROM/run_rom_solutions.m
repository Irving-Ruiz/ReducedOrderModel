clc;
clear;
close all;

% Parameters
Fo = 0:0.01:2;        % Fourier number (dimensionless time)
xi = 1e-5;            % dimensionless position
B  = 10;              % Biot number
a  = inf;             % solution-to-solid volume ratio
m  = 1.5;             % effective dimension
v  = (m - 1)/2;       % Bessel function order parameter
nter = 100;           % number of terms in the series solution
x0 = 1.5;             % initial guess for the first eigenvalue

% Compute eigenvalues
[lambda_n, d_lambda_n] = rom_eigenvalues(nter, B, a, v, x0);

% Compute solutions
Y = rom_average_solution(nter, B, a, v, x0, Fo, lambda_n);
y = rom_local_solution(nter, B, a, v, x0, xi, Fo, lambda_n);

% Plot solutions
figure;

subplot(1,2,1);
semilogy(Fo, Y, 'LineWidth', 1.5);
title('Volume-averaged solution');
xlabel('Fo');
ylabel('\Psi(Fo)');
grid on;

subplot(1,2,2);
semilogy(Fo, y, 'LineWidth', 1.5);
title('Local solution');
xlabel('Fo');
ylabel('\psi(\xi, Fo)');
grid on;