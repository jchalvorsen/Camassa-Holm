function [ U, x, t ] = holdenraynaud(N, T, initial)

%% Configuration
% Domain of x. Note that this is currently locked to [0, 1], but
% hopefully this will be configurable. Worst-case, we transform it to [0,
% 1] and then eventually back.
xmin = 0;
xmax = 1;

%% Preparation
% Spatial step size
h = 1 / N * (xmax - xmin);

% X values in grid
x = xmin + (0:N - 1) * h;

% Find peaks in initial data
peaks = findpeaks(initial(x));

% Determine temporal step size. Use the CFL condition and assume
% the sum of the peaks of the initial data is equal to the velocity of the
% wave (assuming it is a wave).
k = h / abs(sum(peaks));

% t values in grid
t = 0:k:T;

% Amount of time values
M = length(t);

% Determine g
kappa = log( (1 + 2 * N^2 + sqrt(1 + 4*N^2)) / (2 * N ^ 2));
c = 1 / (1 + 2 * N^2 * (1 - exp(-kappa)));
I = 0:N - 1;
g = c * (exp(-kappa * I) + exp(kappa * (I - N))) / (1 - exp(-kappa * N));

% Allocate solution U
U = zeros(M, N);
U(1, :) = initial(x);

% Time solution evaluation run time
ticstart = tic;
progress = 0;
w = waitbar(progress, 'Solving Camassa-Holm...');

%% Execution
for i = 1:M - 1
    u = U(i, :);
    % Extent u on the ends with periodic conditions
    u_periodic = u([end, 1:end, 1]);
    
    % Calculate m
    m = u - fbdiff(u_periodic, h);
    
    % Calculate m_t
    mu = m .* u;
    mu_periodic = mu([end, 1:end, 1]);
    mt = - diff(mu_periodic(1:end - 1)) / h - m .* D(u_periodic, h);
    
    % Applying Euler's Method, we calculate the next "row" of m values
    mnext = m + k * mt;
    
    % Transform mnext back to u
    U(i + 1, :) = ifft(fft(g) .* fft(mnext));
    
    updateProgress( (i + 1) / (M) );
end

elapsed = toc(ticstart);
fprintf('Spent %4.2f seconds on solving equation.\n', elapsed);
close(w);

    function [] = updateProgress(fraction)
        if (fraction > progress * 1.05)
            progress = fraction;
            waitbar(fraction, w);
        end
    end

end

function [ Y ] = fbdiff(X, h)
% Forward, followed by backward finite difference
Y = (X(1:end - 2) - 2 * X(2:end - 1) + X(3:end)) / (h^2);
end

function [ Y ] = D(X, h)
% Average of forward and backward difference
Y =  (X(3:end) - X(1:end - 2))/ (2 * h);
end
