%% Script to test and plot convergence of our implemented method

% Make sure to clear any remnants of U as we'll need the memory
close all
clear all
clear U;

%% Imports
% Make utilities available
path = fileparts(which(mfilename));
addpath(fullfile(path, '../../util'));
clear path;

%% Configuration
% Spatial resolution
I = 6:1:14;
N = 2.^I;
% Maximum time value
T = 15;

xmin = -10;
xmax = 30;

% Compression settings
% nx = number of x values in compressed matrix
nx = 400;
% nt = number of t values in compressed matrix
nt = 400;

%% Initial condition and reference solution
initial = @(x) exp(-abs(x));
%initial = @(x) 2*exp(-abs(x + 5)) + exp(-abs(x - 5));
ref = @(x,t) exp(-abs(x-t));
%ref = @(x,t) 2*exp(-abs(x + 5 -2*t)) + exp(-abs(x - 5 - t));

nM = length(I);
result = zeros(nt,nx,nM);
figure
hold on
color = hsv(nM);
%% Looping to solve equation and check convergence
convergence = zeros(nM,1);
for i = 1:nM
    [U, x, t] = holdenraynaud(N(i), T, [xmin, xmax], initial);
    % generate reference solution:
    reference = zeros(length(t),length(x));
    error = zeros(length(t),1);
    test = zeros(length(t),1);
    for j = 1:length(t)
        reference(j,:) = ref(x, t(j));
        test(j) = norm(reference(j,:),2);
        error(j) = norm(U(j,:) - reference(j,:),2)/norm(reference(j,:),2);
    end
    plot(t,error,'color',color(i,:))
    convergence(i) = error(end) - error(end-1);
    
    M = size(U, 1);
    %% Compression
    % Don't expand matrix if result matrix is smaller
    %nx = min(nx, N(i));
    %nt = min(nt, M);
    [Z, xcomp, tcomp] = compress(x, t, U, nx, nt);
    result(:,:,i) = Z;
end
%% Plotting
%error = result - repmat(reference, [1 1 nM]);
%color = jet(nM);
legenddata = cellstr(num2str((N)'));
ylabel('Relative fault')
xlabel('Time')
legend(legenddata,'Location','NorthWest')

figure
hold on
for i = 1:nM
    plot(xcomp,result(end,:,i)','color',color(i,:))
    %animatedplot(xcomp, tcomp, result(:,:,i));
end
plot(xcomp,ref(xcomp,T),'k')
legenddata(nM + 1 ) = cellstr('Analytical solution');
legend(legenddata,'Location','NorthWest')
xlabel('Travelling peakon')

figure
plot(log(N),log(convergence),'b*-')
hold on
plot(log(N),fliplr(log(N)) - max(log(N)) + max(log(convergence)),'r*-')
legend('Convergence of our method', 'Linear convergence','Location','NorthEast')