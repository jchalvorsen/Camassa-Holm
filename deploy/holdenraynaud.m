%HOLDENRAYNAUD Solves the Camassa-Holm equation by the Holden-Raynaud scheme.
%   [U, x, t] = HOLDENRAYNAUD(N, T, xdomain, initial) solves the 
%   Camassa-holm equation on the spacial domain xdomain = [xmin, xmax] and 
%   temporal interval [0, T] with a resolution of N points on the x-axis, 
%   using the function initial(x) as a periodic initial condition. x and t 
%   are row vectors of spatial and temporal values respectively, with 
%   lengths N = length(x) and M = length(t), such that the result matrix U
%   has dimensions M x N. Note that the scheme used here is extended with
%   the proposition put forward by M. Dahlby in order to handle antipeakons.
% 
%   HOLDENRAYNAUD also accepts additional optional arguments. Following
%   standard conventions, the options are supplied in a (..., 'key', value)
%   manner. Available options:
%       - 'ShowProgress': true/false. Toggles whether or not to show a
%                         progressbar as the calculations are running.
%                         Default: true.
%       - 'PrintTiming':  true/false. Toggles whether or not to print the
%                         computational time when finished. Default: true.
%       - 'dt':           double. Override the time step given from the CFL
%                         condition heuristic with a user-supplied value.
%                         Not recommended for ordinary use.
%  Examples:
%
%  initial = @(x) exp(-abs(x));
%  [U, x, t] = holdenraynaud(1024, 10, [-10, 10], initial, 'ShowProgress', false);
%  surf(x, t, U);

function [ U, x, t ] = holdenraynaud(N, T, xdomain, initial, varargin)

[ showprogress, printtiming, timestep ] = parse(varargin);

%% Sanity checks
assert(isvector(xdomain), 'xdomain must be a vector.');
assert(length(xdomain) == 2, 'xdomain must have length 2.');
assert(xdomain(2) > xdomain(1), 'xdomain = [XMIN, XMAX] must have XMAX > XMIN');
assert(N > 1, 'N must be greater than 1.');
assert(timestep >= 0, 'Supplied timestep must be non-negative');

%% Configuration
% Domain of x.
xmin = xdomain(1);
xmax = xdomain(2);

%% Preparation
% Spatial step size
h = (xmax - xmin) / N;

% x values in grid
x = xmin + (0:N-1) * h;

% Generate initial data
initialdata = initial(x);

% Temporal step size
dt = calculatetimestep(initial, xmin, xmax, h);
if timestep == 0
    % Use time step calculcated from CFL
    k = dt;
else
    % Use user-supplied time step
    k = timestep;
    if k > dt
       warning(['User-supplied time step is larger than CFL-based heuristic. ' ...
           'Stability may suffer as a result.']);
    end
end

% Adjust time step (downward) so that t(end) == T. M is number of time values
M = fix(T/ k) + 2;
k = T / (M - 1);

% t values in grid
t = 0:k:T;

% Determine g (see paper by Holden, Raynaud), and pre-calculate fft(g).
% Note that by replacing n in their derivation by K = 1 / h, 
% the method works on arbitrary intervals.
K = 1 / h;
kappa = log( (1 + 2 * K^2 + sqrt(1 + 4*K^2)) / (2 * K ^ 2));
c = 1 / (1 + 2 * K^2 * (1 - exp(-kappa)));
I = 0:N - 1;
g = c * (exp(-kappa * I) + exp(kappa * (I - N))) / (1 - exp(-kappa * N));
fftg = transpose(fft(g));

% Allocate solution U
U = zeros(M, N);
U(1, :) = initialdata;

% Start timer and show progressbar if enabled
ticstart = tic;
progress = 0;
if showprogress
    w = waitbar(progress, 'Solving Camassa-Holm...');
end

% Difference operator matrices
B = backward_matrix(N, h);
F = forward_matrix(N, h);
C = (backward_matrix(N, h) + forward_matrix(N, h)) / 2;
A = speye(N) - B * F;

%% Execution
for i = 1:M - 1
    % Calculate m_t by the modified Holden-Raynaud scheme
    u = transpose(U(i, :));
    m = A * u;
    mt = - (B * (m .* max(u, 0)) + F * (m .* min(u, 0))) - m .* (C * u);
    
    % Calculate new values for m using Euler's method
    m = m + k * mt;
    
    % Transform m back to u by convolution
    U(i + 1, :) = ifft(fftg .* fft(m));
    updateProgress( (i + 1) / (M) );
end

%% Cleanup and return
elapsed = toc(ticstart);
if printtiming
    fprintf('Solved Camassa-Holm in %4.2f seconds with a spatial grid size of %d.\n', ...
        elapsed, N);
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
function [ B ] = backward_matrix(N, h)
B = spdiags([- ones(N, 1), ones(N, 1) ], -1:0, N, N) / h;
B(1,N) = -1 / h;
end

function [ F ] = forward_matrix(N, h)
F = - transpose(backward_matrix(N, h));
end

%% Timestep and CFL condition
function [ k ] = calculatetimestep(initial, a, b, h)
% Calculate potential maximum area of a peakon to estimate maximum height
areaAbove = integral(@(x) max(initial(x), 0), a, b);
areaBelow = - integral(@(x) min(initial(x), 0), a, b);
A = max(areaAbove, areaBelow);
c = A / 2;

% Use the CFL condition and assume the maximum height is equal to the
% velocity of the wave (assuming it is a wave).
% Multiply by a heuristic "goodness factor" to overestimate its maximum velocity.
k = h / (1.1 * c);

end

%% Parameter parsing for holdenraynaud
function [ showprogress, printtiming, timestep ] = parse(options)
% Parses additional options to holdenraynaud

% Set default values for options
showprogress = true;
printtiming = true;
timestep = 0;

count = length(options);
for k = 1:2:count
    % Sanity checks
    parameter = lower(options{k});
    missingMessage = strcat('Missing parameter value for parameter ''', ...
        parameter, '''.');
    assert(k + 1 <= count, missingMessage);
    
    value = options{k + 1};
    assert(~isempty(value), missingMessage);
    
    % Note: Lower-case for case insensitivity
    switch parameter
        case 'showprogress'
            assertlogical(parameter, value);
            showprogress = value;
        case 'printtiming'
            assertlogical(parameter, value);
            printtiming = value;
        case 'dt'
            assertnumber(parameter, value);
            timestep = value;
    end
    
    
end

    function assertlogical(parameter, value)
        assert(islogical(value), strcat('Parameter value for parameter ''', ...
            parameter, ''' is not logical (true/false).'));
    end

    function assertnumber(parameter, value)
        assert(isscalar(value) && isnumeric(value), strcat('Parameter value for parameter ''', ...
            parameter, ''' is not a numerical scalar.'));
    end
end
