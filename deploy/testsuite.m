%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Testsuite for solving the Camassa-Holm equation 
% by the (modified) Holden-Raynaud scheme
% -------------------------------------------------
% 
% © 2014 J. C. Halvorsen, T. Østhus, A. Longva.
%   
% The following script should be fairly self-explanatory, but there are
% particularly two utility functions we have implemented that may be worth
% elaborating on:
% 
% - compress: This function compresses (or expands) an arbitrarily-sized
%             matrix into an arbitrarily configured size by interpolation. 
%             Very convenient when plotting high-resolution datasets, 
%             particularly for 3D plots, as it retains the ability to 
%             rotate and manipulate the viewport smoothly, while still 
%             maintaining the structure of the high-resolution matrix.
%
% - animatedplot: A simple function that animates 2D plots, taking an 
%             arbitrary number of T x X matrices as input, each 
%             representing a data series. 
%
% We have provided a selection of initial conditions that show interesting
% effects. You're also welcome to create your own, but note that our CFL
% condition is not sufficient for stability, only necessary, meaning you
% may be able to find certain classes of initial data for which the scheme
% is not stable.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

% Compression settings (for plotting)
% nx = number of x values in compressed matrix
nx = 400;
% nt = number of t values in compressed matrix
nt = 300;

%% Initial condition (as a function of x)
% Define a basic peakon for use in initial functions
peakon = @(x, t) exp(-abs(x - t));

% Comment/uncomment to pick your desired initial condition. Note that
% initial conditions should be periodic or approximately periodic.

% Single peakon
xmin = -20;
xmax = 20;
initial = @(x) peakon(x, 0);

% Double peakon
% xmin = -20;
% xmax = 20;
% initial = @(x) peakon(x, -5) + peakon(x, 5);

% Peakon-peakon interaction (may require very fine resolution)
% xmin = -20;
% xmax = 20;
% initial = @(x) peakon(x, -5) + 0.7 * peakon(x, 3);

% Peakon-antipeakon interaction
% xmin = -20;
% xmax = 20;
% initial = @(x) peakon(x, -5) - peakon(x, 5);

% Smooth, periodic initial condition
% xmin = -20;
% xmax = 20;
% initial = @(x) sin(2 * pi * x / (xmax - xmin)) + 1;

% Smooth initial condition (not periodic, but interesting nonetheless)
% xmin = -20;
% xmax = 20;
% initial = @(x) sin(pi * x / (xmax - xmin)) + 1;

% Highly oscillatory smooth initial condition
% xmin = 0;
% xmax = 40;
% initial = @(x) sin(pi * x) + 1;

% Smooth soliton solution shape of the KdV equation
% xmin = -20;
% xmax = 20;
% initial = @(x) sech(x / 2).^2 / 2;

% Custom initial function
% xmin = ???
% xmax = ???
% initial = @(x) ????'

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
xlabel('x')
ylabel('t')
animatedplot(xcomp, tcomp, Z)