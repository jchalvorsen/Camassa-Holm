% Make sure to clear any remnants of U as we6'll need the memory
clear U;

%% Imports
% Make utilities available
path = fileparts(which(mfilename));
addpath(fullfile(path, '../../util'));
clear path;

%% Configuration
% Spatial resolution
N = 1024;
% Maximum time value
T = 240;

xmin = -10;
xmax = 60;

% Compression settings
% nx = number of x values in compressed matrix
nx = 500;
% nt = number of t values in compressed matrix
nt = 600;

%% Initial condition (as a function of x)
a = 1;
initial = @(x) exp(-abs(x));

%% Solve equation
[U, x, t, V] = linearholdenraynaud(N, T, [xmin, xmax]);
M = size(U, 1);

%% Compression
% Don't expand matrix if result matrix is smaller
nx = min(nx, N);
nt = min(nt, M);
[Z, xcomp, tcomp] = compress(x, t, U, nx, nt);
V = compress(x, t, V, nx, nt);

%% Plotting
figure
surf(xcomp, tcomp, Z)
shading flat;
animatedplot(xcomp, tcomp, Z, V)
xlabel('x')
ylabel('y')