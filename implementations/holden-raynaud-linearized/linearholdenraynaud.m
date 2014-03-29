function [ U, x, t ] = holdenraynaud(N, T, xdomain, initial, varargin)

[ showprogress, printtiming ] = parse(varargin);

%% Sanity checks
assert(isvector(xdomain), 'xdomain must be a vector.');
assert(length(xdomain) == 2, 'xdomain must have length 2.');
assert(xdomain(2) > xdomain(1), 'xdomain = [XMIN, XMAX] must have XMAX > XMIN');
assert(N > 1, 'N must be greater than 1.');

%% Configuration
% Domain of x. 
xmin = xdomain(1);
xmax = xdomain(2);

%% Linearization constants
u0 = @(y) cosh(min(y, 1 - y));

%% Preparation
% Spatial step size
h = 1 / N;

% y and corresponding x values in grid
y = (0:N - 1) * h;
x = (xmax - xmin) * y + xmin;

% Generate initial data
initialdata = initial(x);

% Find peaks in initial data. Add boundaries if necessary
peaks = findpeaks(initialdata);
if initialdata(end) > initialdata(end - 1); peaks = [ peaks initialdata(end) ]; end
if initialdata(1) > initialdata(2); peaks = [ initialdata(1) peaks ]; end

% Determine temporal step size. Use the CFL condition and assume
% the sum of the peaks of the initial data is equal to the velocity of the
% wave (assuming it is a wave).
k = h / abs(sum(peaks));

% t values in grid
t = 0:k:T;

% Amount of time values
M = length(t);

% Determine g (see paper by Holden, Raynaud), and pre-calculate fft(g)
kappa = log( (1 + 2 * N^2 + sqrt(1 + 4*N^2)) / (2 * N ^ 2));
c = 1 / (1 + 2 * N^2 * (1 - exp(-kappa)));
I = 0:N - 1;
g = c * (exp(-kappa * I) + exp(kappa * (I - N))) / (1 - exp(-kappa * N));
fftg = fft(g);

% Note: While the x domain can be any sensible domain on the real line,
% we only work with the interval [0, 1]. Hence, we transform x to the
% variable y, which resides on the unit interval. We accomplish this simply
% by letting the initial data for the x values on the interval [xmin, xmax]
% be equal to the initial data for y on [0, 1]

% Allocate solution U
U = zeros(M, N);
U(1, :) = initialdata;

% Time solution evaluation run time
ticstart = tic;
progress = 0;
if showprogress
    w = waitbar(progress, 'Solving Camassa-Holm...');
end

%% Linearization
U0 = u0(y);
M0 = U0 - fbdiff(U0([end, 1:end, 1]), h);

A0 = bdiff(M0, h);
B0 = (3 / 2) * M0;
C0 = M0 / 2;
D0 = (3 / 2) * bdiff(U0, h) + (1 / 2) * fdiff(U0, h);
E0 = U0;

alpha = D0 .* M0 + E0 .* A0;

%% Execution
for i = 1:M - 1
    u = U(i, :);
    % Extend u on the ends with periodic conditions
    u_periodic = u([end, 1:end, 1]);
    
    % Calculate m
    m = u - fbdiff(u_periodic, h);
    
    % Calculate m_t
    %mu = m .* u;
    %mu_periodic = mu([end, 1:end, 1]);
    %mt = - diff(mu_periodic(1:end - 1)) / h - m .* D(u_periodic, h);
    mt = alpha - A0 .* u - B0 .* bdiff(u, h) - C0 .* fdiff(u, h) - D0 .* m - E0 .* bdiff(m, h);
    
    % Applying Euler's Method, we calculate the next "row" of m values
    mnext = m + k * mt;
    
    % Transform mnext back to u by convolution
    U(i + 1, :) = ifft(fftg .* fft(mnext));
    
    updateProgress( (i + 1) / (M) );
end

%% Exit

elapsed = toc(ticstart);
if printtiming
    fprintf('Spent %4.2f seconds on solving equation.\n', elapsed);
end

if showprogress
    close(w);
end

    function [] = updateProgress(fraction)
        if showprogress
            if (fraction - progress >= 0.05 || progress >= 1.0)
                progress = fraction;
                waitbar(fraction, w);
            end
        end
    end
end

%% Finite differences
function [ Y ] = fdiff(X, h)
Y = diff(X([1:end, 1])) / h;
end

function [ Y ] = bdiff(X, h)
Y = diff(X([end, 1:end])) / h;
end

function [ Y ] = fbdiff(X, h)
% Forward, followed by backward finite difference
Y = (X(1:end - 2) - 2 * X(2:end - 1) + X(3:end)) / (h^2);
end

function [ Y ] = D(X, h)
% Average of forward and backward difference
Y =  (X(3:end) - X(1:end - 2))/ (2 * h);
end

%% Parameter parsing for holdenraynaud
function [ showprogress, printtiming ] = parse(options)
% Parses additional options to holdenraynaud

% Set default values for options
showprogress = true;
printtiming = true;

count = length(options);
for k = 1:2:count
    % Sanity checks
    parameter = lower(options{k});
    missingMessage = strcat('Missing parameter value for parameter ''', ...
        parameter, '''.');
    assert(k + 1 <= count, missingMessage);
    
    value = options{k + 1};
    assert(~isempty(value), missingMessage);
    assert(islogical(value), strcat('Parameter value for parameter ''', ...
        parameter, ''' is not logical (true/false).'));
    
    % Note: Lower-case for case insensitivity
    switch parameter
        case 'showprogress'
            showprogress = value;
        case 'printtiming'
            printtiming = value;
    end
end
end