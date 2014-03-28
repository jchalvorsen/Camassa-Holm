% Make sure to clear any remnants of U as we'll need the memory
clear U;

%% Imports
% Make utilities available
path = fileparts(which(mfilename));
addpath(fullfile(path, '../../util'));
clear path;

%% Configuration
% Spatial resolution
N = 2048;
% Maximum time value
T = 5;

xmin = -10;
xmax = 20;

% Compression settings
% nx = number of x values in compressed matrix
nx = 100;
% nt = number of t values in compressed matrix
nt = 600;

%% Initial condition
initial = @(x) 2*exp(-abs(x));

%% Solve equation
[U, x, t] = holden_matrix_method(N, T, [xmin, xmax], initial);
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