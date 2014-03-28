function [] = forwardtime_centralspace()

%% Configuration
M = 128;
xmin = -30;
xmax = 30;
tmax = 400;

% Initial condition
g = @(x) exp(-abs(x-2)) + exp(-abs(x-4));
g = @(x) exp(-(x-2).^2);
%g = @(x) 10 ./ ((3 + abs(x)) .^ 2);
%g = @(x) abs(x) ./ max(x);

%% Preparation
h = (xmax - xmin) / (M - 1)
k = h^2%min(tmax / 100, h^2)

X = xmin:h:xmax;
T = [0, k:k:tmax ];

%% Execution

U = zeros(length(T), length(X));
U(1, :) = g(X);

% Calculate U(n + 1, :) in every loop iteration
for n = 1:length(T) - 1  
    row = U(n, :);
    
    % Calculate the difference terms of the right hand side
    d = central([row(end-1), row, row(2)], h);
    d2 = central2([row(end-1), row, row(2)], h);
    d3 = central3([row(end-2:end-1), row, row(2:3)], h);

    % Solve the system of equations
    A = build_matrix(length(X), h);
    F = transpose(row - d2 + (k * d2 .* (2 * d)) + k * row .* (d3 - 3 * d));
    U(n + 1, :) = A \ F;
    
    % Scale
%     U(n + 1, :) = (U(n + 1, :) - min(U(n + 1, :)));
%     U(n + 1, :) = U(n + 1, :) / max(U(n + 1, :));
    
end

figure
plot(X, g(X))

animated_plot(X, U, 0.05)
 
% figure;
% surf(X, T, U)
% shading flat;
% xlabel('x')
% ylabel('t')

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
ylim([min(U(1, :)), 2 * max(U(1, :))]);
pause(step)
for i = 2:size(U, 1)
    y = U(i, :);
    set(h,'YData', y)
    drawnow
    pause(step);
end

end