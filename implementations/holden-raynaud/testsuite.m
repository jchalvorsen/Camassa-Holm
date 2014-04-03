% Make sure to clear any remnants of U as we'll need the memory
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
T = 50;
% Domain
xmin = -20;
xmax = 20;

% Compression settings (for plotting)
% nx = number of x values in compressed matrix
nx = 400;
% nt = number of t values in compressed matrix
nt = 200;

%% Initial condition (as a function of x)
peakon = @(x, t) exp(-abs(x - t));

% Comment/uncomment to pick your desired initial condition. Note that
% initial conditions must be periodic or approximately periodic.

% Single peakon
%initial = @(x) peakon(x, 0);

% Double peakon
%initial = @(x) peakon(x, -5) + peakon(x, 5);

% Peakon-peakon interaction (may require very fine resolution)
%initial = @(x) peakon(x, -5) + 0.7 * peakon(x, 3);

% Peakon-antipeakon interaction
%initial = @(x) peakon(x, -5) - peakon(x, 5);

% Smooth initial condition
initial = @(x) sin(pi * x / (xmax - xmin)) + 1;

% Custom initial function
%initial = @(x) ????
%% Solve equation
[U, x, t] = holdenraynaud(N, T, [xmin, xmax], initial, 'showprogress', true);
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