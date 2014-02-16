function [] = forwardtime_centralspace()

%% Configuration
M = 256;
xmin = 0;
xmax = 4;
tmax = 0.1;

% Initial condition
g = @(x) exp(-abs(x-2));

%% Preparation
h = (xmax - xmin) / (M - 1)
k = h^2

X = xmin:h:xmax;
T = [0, k:k:tmax ];

%% Execution

U = zeros(length(T), length(X));
U(1, :) = g(X);

% Calculate U(n + 1, :) in every loop iteration
for n = 1:length(T) - 1  
    row = U(n, :);
    
    % Calculate the difference terms of the right hand side
    d = central([row(end), row, row(1)], h);
    d2 = central2([row(end), row, row(1)], h);
    d3 = central3([row(end-1:end), row, row(1:2)], h);

    % Solve the system of equations
    A = build_matrix(M, h);
    F = transpose(d2 .* (1 + 2 * d) + row .* (1 + k * d3 - 3 * k * d));
    U(n + 1, :) = A \ F;
end

plot(X, g(X))

figure
animated_plot(X, U, 0.1)

% figure;
% surf(X, T, U)
% xlabel('x')
% ylabel('t')
% shading flat;

end

function [ T ] = central(V, d)
T = (V(3:end) - V(1:end-2)) / (2 * d);
end

function [ T ] = central2(V, d)
T = (V(3:end) - 2 * V(2:end-1) + V(1:end-2)) / (d ^ 2);
end

function [ T ] = central3(V, d)
T = (V(5:end) - 2 * V(4:end-1) + 2 * V(2:end-3) - V(1:end-4)) / (2 * d^3);
end

function [ A ] = build_matrix(M, h)
neighbours = (- 1 / h^2) * ones(M, 1);
main = (2 / h^2 + 1) * ones(M, 1);
A = spdiags([neighbours, main, neighbours], -1:1, M, M);
A(1, M) = (- 1 / h^2);
A(M, 1) = (- 1 / h^2);
end

function animated_plot(x, U, step)
figure
y = U(1, :);
h = plot(x, y);
xlim([min(x), max(x)])
ylim([max(min(U)), min(max(U))])
pause(step)
for i = 2:size(U, 1)
    y = U(i, :);
    set(h,'YData', y)
    drawnow
    pause(step);
end

end