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

%% Preparation
% Spatial step size
h = (xmax - xmin) / N;
%h = 1 / N;

% y and corresponding x values in grid
%y = (0:N - 1) * h;
%x = (xmax - xmin) * y + xmin;
x = xmin + (0:N-1) * h;

% Generate initial data
initialdata = initial(x);

% Find peaks in initial data. Add boundaries if necessary
peaks = findpeaks(initialdata);
if initialdata(end) > initialdata(end - 1); peaks = [ peaks initialdata(end) ]; end
if initialdata(1) > initialdata(2); peaks = [ initialdata(1) peaks ]; end

% Determine temporal step size. Use the CFL condition and assume
% the sum of the peaks of the initial data is equal to the velocity of the
% wave (assuming it is a wave). Multiply by a factor to overestimate its
% maximum velocity
k = h / (1.05 * abs(sum(peaks)));

% t values in grid
t = 0:k:T;

% Amount of time values
M = length(t);

K = 1 / h;
% Determine g (see paper by Holden, Raynaud), and pre-calculate fft(g)
kappa = log( (1 + 2 * K^2 + sqrt(1 + 4*K^2)) / (2 * K ^ 2));
c = 1 / (1 + 2 * K^2 * (1 - exp(-kappa)));
I = 0:N - 1;
g = c * (exp(-kappa * I) + exp(kappa * (I - N))) / (1 - exp(-kappa * N));
fftg = transpose(fft(g));

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

% Difference operator matrices
% Forward + backward?
e = ones(N,1);
A = spdiags([-e/h^2, (1 + 2/h^2) * e, -e/h^2], -1:1, N, N);
A(N,1) = -1/h^2;
A(1,N) = -1/h^2;
% B - backward difference
B = spdiags([-e/h, e/h], -1:0, N, N);
B(1,N) = -1/h;

% F - forward difference
F = -B';

% C - central difference
C = spdiags([-e/(2*h), e/(2*h)], [-1, 1], N, N);
C(1,N) = -1/(2*h);
C(N,1) = 1/(2*h);

%% Execution
for i = 1:M - 1
    u = transpose(U(i, :));
    % Extend u on the ends with periodic conditions
    u_periodic = u([end, 1:end, 1]);
    
    % Calculate m
    m = u - fbdiff(u_periodic, h);
    %m = A * u;
    
    % Calculate m_t
    mu = m .* u;
    %mt = - bdiff(mu, h) - m .* D(u_periodic, h);
    mt = - (B*(m.*max(u,0)) + F*(m.*min(u,0))) - m.*(C*u);
    
    % Applying Euler's Method, we calculate the next "row" of m values
    mnext = m + k * mt;
    
    % Transform mnext back to u by convolution
    U(i + 1, :) = ifft(fftg .* fft(mnext));
    %U(i + 1, :) = (A \ mnext)';
    
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
Y = (X(3:end) - X(1:end - 2))/ (2 * h);
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
