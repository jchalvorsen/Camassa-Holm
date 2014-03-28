function [ U, x, t ] = holden_matrix_method(N, T, xdomain, initial, varargin)

[ showprogress, printtiming ] = parse(varargin);

%% Sanity checks
assert(isvector(xdomain), 'xdomain must be a vector.');
assert(length(xdomain) == 2, 'xdomain must have length 2.');
assert(xdomain(2) > xdomain(1), 'xdomain = [XMIN, XMAX] must have XMAX > XMIN');
assert(N > 1, 'N must be greater than 1.');

%% Configuration
% Domain of x. Note that this is currently locked to [0, 1], but
% hopefully this will be configurable. Worst-case, we transform it to [0,
% 1] and then eventually back.
xmin = xdomain(1);
xmax = xdomain(2);

%% Preparation
% Spatial step size
h = 1 / N * (xmax - xmin);

% X values in grid
x = xmin + (0:N - 1) * h;

% Generate initial data
initialdata = initial(x);

% Find peaks in initial data. Add boundaries if necessary
peaks = findpeaks(initialdata);
if initialdata(end) > initialdata(end - 1); peaks = [ peaks initialdata(end) ]; end
if initialdata(1) > initialdata(2); peaks = [ initialdata(1) peaks ]; end

% Determine temporal step size. Use the CFL condition and assume
% the sum of the peaks of the initial data is equal to the velocity of the
% wave (assuming it is a wave).
k = h / (2*abs(sum(peaks)));

% t values in grid
t = 0:k:T;

% Amount of time values
M = length(t);

% Allocate solution U
U = zeros(M, N);
U(1, :) = initialdata;

% Time solution evaluation run time
ticstart = tic;
progress = 0;
if showprogress
    w = waitbar(progress, 'Solving Camassa-Holm...');
end

%% Execution
% A - operator to go from u to m
e = ones(N,1);
A = spdiags([-e/h^2, (1 + 2/h^2) * e, -e/h^2], -1:1, N, N);
A(N,1) = -1/h^2;
A(1,N) = -1/h^2;

% B - backwards difference
B = spdiags([-e/h, e/h], -1:0, N, N);
B(1,N) = -1/h;

% C - central difference
C = spdiags([-e/(2*h), e/(2*h)], [-1, 1], N, N);
C(1,N) = -1/(2*h);
C(N,1) = 1/(2*h);


for i = 1:M-1
    u = U(i,:)';
    % updating waitbar
    if showprogress
        updateProgress( (i + 1) / (M) );
    end
    % TEST OF NEW FORMAT
    m = A*u;
    mt = m + k*(- B*(m.*u) - m.*(C*u) );
    U(i+1,:) = (A\mt)';
    % Everythin in one operation below
    %u = u + k * (A\((-B*((A*u).*u)) - (A*u).*(C*u)));
end

%% End execution

elapsed = toc(ticstart);
if printtiming
    fprintf('Spent %4.2f seconds on solving equation.\n', elapsed);
end

if showprogress
    close(w);
end

    function [] = updateProgress(fraction)
        if (fraction - progress >= 0.05 || progress >= 1.0)
            progress = fraction;
            waitbar(fraction, w);
        end
    end
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
