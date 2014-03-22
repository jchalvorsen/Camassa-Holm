function [ U ] = holdenraynaud()

%% Configuration
% Spatial resolution
N = 512;
% Domain of x
xmin = 0;
xmax = 1;
% Maximum time value
T = 1;
% Initial condition (as a function of x)
a = 1 / (xmax - xmin);  
initial = @(x) cosh(min(x, a - x)) / sinh(a / 2);

% Compression settings
% nx = number of x values in compressed matrix
nx = 80;
% nt = number of t values in compressed matrix
nt = 500;

%% Preparation
% Spatial step size
h = 1 / N * (xmax - xmin);

% X values in grid
x = (0:N - 1) * h;

% Determine temporal step size. Use the CFL condition and assume
% the maxmimum size of the initial data is equal to the velocity of the
% wave (assuming it is a wave).
k = h / max(abs(initial(x)));

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
w = waitbar(0, 'Solving Camassa-Holm...');

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
close(w)

%% Compression

ticstart = tic;
[Z, xcomp, tcomp] = compress(x, t, U, nx, nt);
compresselapsed = toc(ticstart);

fprintf('Spent %4.2f seconds on compressing solution before plotting.\n', ...
    compresselapsed);

%% Plotting

figure
%surf(xcomp, tcomp, Z)
shading flat;
animatedplot(xcomp, tcomp, Z)
xlabel('x')
ylabel('y')

    function [] = updateProgress(fraction)
        waitbar(fraction, w);
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

function [ ] = animatedplot(x, t, U)
% Target framerate
FRAMERATE = 30;
FRAMELENGTH = 1 / FRAMERATE;

% Plot initial data
figure
plothandle = plot(x, U(1, :));
Umax = max(U(1, :));
Umin = min(min(U));
ylim([Umin, Umax]);
for j = 2:length(t)
    set(plothandle, 'YData', U(j, :))
    drawnow
    pause(FRAMELENGTH)
end

end