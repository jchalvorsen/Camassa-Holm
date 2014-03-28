% Make sure to clear any remnants of U as we'll need the memory
clear U;

%% Imports
% Make utilities available
path = fileparts(which(mfilename));
addpath(fullfile(path, '../../util'));
clear path;

%% Configuration
% Spatial resolution
N = 256;
% Maximum time value
T = 5;

xmin = 0;
xmax = 1;

% Compression settings
% nx = number of x values in compressed matrix
nx = 100;
% nt = number of t values in compressed matrix
nt = 600;

%% Initial condition (as a function of x)
a = 1;
initial = @(x) cosh(min(x, a - x));
%initial = @(x) [cosh(min(x(1:floor(end/2)), a - x(1:floor(end/2)))), ...
%    repmat(cosh(min(x(floor(end / 2 + 1)), a - x(floor(end / 2 + 1)))), 1, ...
%    length(x) - floor(length(x) / 2))] - ...
%    cosh(min(x(floor(end / 2 + 1)), a - x(floor(end / 2 + 1))));
%initial = @(x) cosh(min(x, a - x)) + circshift(0.5 * cosh(min(x, a - x)), ...
%    repmat(round(length(x) / 2), length(x), 1));
%initial = @(x) exp(-10 * abs(x - 1.2)) + exp(-10 * abs(x - 0.2));

%% Solve equation
[U, x, t] = linearholdenraynaud(N, T, [xmin, xmax], initial);
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