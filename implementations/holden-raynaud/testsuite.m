% Make sure to clear any remnants of U as we'll need the memory
clear U;

%% Imports
% Make utilities available
path = fileparts(which(mfilename));
addpath(fullfile(path, '../../util'));
clear path;

%% Configuration
% Spatial resolution
N = 5000;
% Maximum time value
T = 5;

xmin = -10;
xmax = 10;

% Compression settings
% nx = number of x values in compressed matrix
nx = 500;
% nt = number of t values in compressed matrix
nt = 600;

%% Initial condition (as a function of x)
a = 1;
%initial = @(x) cosh(min(x, a - x)) / sinh(a / 2);

%initial = @(x) [cosh(min(x(1:floor(end/2)), a - x(1:floor(end/2)))), ...
%    repmat(cosh(min(x(floor(end / 2 + 1)), a - x(floor(end / 2 + 1)))), 1, ...
%    length(x) - floor(length(x) / 2))] - ...
%    cosh(min(x(floor(end / 2 + 1)), a - x(floor(end / 2 + 1))));
%initial = @(x) cosh(min(x, a - x)) + circshift(0.5 * cosh(min(x, a - x)), ...
%    repmat(round(length(x) / 2), length(x), 1));
%initial = @(x) exp(-10 * abs(x - 1.2)) + exp(-10 * abs(x - 0.2));
%initial = @(x) cosh(min(x, a - x)) / sinh(a / 2);
%initial = @(x) 3 * exp(-abs(x)) - 3 * exp(-abs(x - 10));
%initial = @(x) exp(-abs(x + 2)) + exp(-abs(x - 2));
%initial = @(x) 5 * exp(-abs(x + 5)) + 5 * exp(-abs(x - 25));
%initial = @(x) exp(-abs(x + 5)) - exp(-abs(x - 5));
%initial = @(x) sin(2 * pi * x / 20) + 1;
%initial = @(x) sin(pi * x / 2);
initial = @(x) exp(-abs(x));
%initial = @(x) sin(pi * x / 10);
% m1 = (exp(-5) + 6) / (2 * exp(-5) + 3);
% m2 = (exp(-5) + 2/3) / (exp(-5) + 3);
% x1 = log(18 * exp(-10) / (exp(-5)) + 6);
% x2 = log(40 * exp(-10) + 60 * exp(-5));
% initial = @(x) m1 * exp(-abs(x - x1)) + ...
%      m2 * exp(-abs(x - x2));

%% Solve equation
[U, x, t] = holdenraynaud(N, T, [xmin, xmax], initial);
M = size(U, 1);

%% Compression
% Don't expand matrix if result matrix is smaller
nx = min(nx, N);
nt = min(nt, M);
[Z, xcomp, tcomp] = compress(x, t, U, nx, nt);

%% Plotting
figure
surf(xcomp, tcomp, Z)
shading flat;
animatedplot(xcomp, tcomp, Z)
xlabel('x')
ylabel('y')