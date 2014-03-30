% Make sure to clear any remnants of U as we'll need the memory
close all
%clear all
clear U;

%% Imports
% Make utilities available
path = fileparts(which(mfilename));
addpath(fullfile(path, '../../util'));
clear path;

%% Configuration
% Spatial resolution
I = 8:1:10;
N = 2.^I;
% Maximum time value
T = 10;

xmin = -10;
xmax = 20;

% Compression settings
% nx = number of x values in compressed matrix
nx = 100;
% nt = number of t values in compressed matrix
nt = 100;

%% Initial condition and building reference solution
initial = @(x)  exp(-abs(x));% - 2*exp(-abs(x - 5));% + exp(-abs(x - 5));
ref = @(x,t) 2*exp(-abs(x-2*t));
reference = zeros(nt,nx);
Ts = 0 : T/nt : T-T/nt;
for i = 1:nt
    reference(i,:) = ref(xmin:(xmax-xmin)/nx:xmax-1/nx, Ts(i));
end

nM = length(I);
result = zeros(nt,nx,nM);
%% Looping to solve equation and check convergence
for i = 1:nM
    [U, x, t] = holden_matrix_linear(N(i), T, [xmin, xmax], initial);
    M = size(U, 1);
    
    %% Compression
    % Don't expand matrix if result matrix is smaller
    nx = min(nx, N(i));
    nt = min(nt, M);
    [Z, xcomp, tcomp] = compress(x, t, U, nx, nt);
    result(:,:,i) = Z;
end
%% Plotting
error = result - repmat(reference, [1 1 nM]);
%color = jet(nM);

for i = 1:nM
    %norms(i) = norm(error(end,:,i)-error(end-1,:,i))/norm(reference(end,:));
    
    animatedplot(xcomp, tcomp, result(:,:,i));
end
%plot(log(N),log(norms),'*-r')
%hold on
%plot(log(N),log(fliplr(N.^(1/2))), '*-b')