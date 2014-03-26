clear all;

%% Imports
% Make utilities available
path = fileparts(which(mfilename));
addpath(fullfile(path, '../../util'));
clear path;

%% Purpose of this script
% Generating very high-resolution solutions requires huge amounts of
% memory. To get around this issue, we solve our equation in a stepwise
% process, "compressing" the results for each step, and using the end
% results of the previous iteration as the initial conditions for the next.
% 
% By "compression", we mean linearly interpolating a matrix which
% represents a grid of numerical data in such a way that the "structure" of
% the data is approximately intact (very accurately for high-resolution
% data), while the size of the compressed matrix is dramatically reduced,
% and thus suitable for interactive plotting.

%% Configuration
% Spatial resolution (number of points along the x-axis)
N = 1000;

% Solve for the total duration [0, T]
T = 10;

% Stepsize
S = 0.1;

% Compression settings
nx = 200;
nt = 1000;

%% Initial condition (as a function of x)
a = 1;
%initial = @(x) cosh(min(x, a - x));
initial = @(x) cosh(min(x, a - x)) + circshift(0.5 * cosh(min(x, a - x)), repmat(round(length(x) / 3), length(x), 1));

%% Preparation
Z = zeros(0, 0);
nsteps = ceil(T / S);

cx = [];
ct = [];

%% Execution
for i = 1:nsteps
   [U, x, t] = holdenraynaud(N, S, initial);
   [V, vx, vt] = compress(x, t, U, nx, nt);
   Z = [Z; V];
   
   cx = vx;
   
   if ~isempty(ct)
        ct = [ct, (ct(end) + vt)];
   else
        ct = vt;
   end
   
   % Redefine initial to use last value in U
   initial = @(x) U(end, :);
end

animatedplot(cx, ct, Z);

