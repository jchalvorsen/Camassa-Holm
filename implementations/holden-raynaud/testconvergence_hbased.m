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
I = 5:16;
h = zeros(1, length(I));

% Maximum time value
T = 5;
xmin = -10;
xmax = 15;

% Compression settings
% nx = number of x values in compressed matrix
nx = 400;
% nt = number of t values in compressed matrix
nt = 400;

%% Initial condition and reference solution
ref = @(x,t) exp(-abs(x-t));
initial = @(x) ref(x, 0);

% Spatial steps

% result = zeros(nt,nx,nM);
% figure
% hold on
% color = hsv(nM);
%% Loop to solve equation and check convergence
errors = zeros(1, length(I)); 
for j = 1:length(I);
    N = 2 ^ I(j);
    [U, x, t] = holdenraynaud(N, T, [xmin, xmax], initial);
    
    h(j) = x(2) - x(1);
    % generate reference solution:
    reference = zeros(length(t),length(x));
    for l = 1:length(t)
        reference(l,:) = ref(x, t(l));
    end
    
    % Error function norm
    errors(j) = sqrt(h(j)) * norm(U(l,:) - reference(l,:),2);
    
    %plot(t,error,'color',color(i,:))
    
    %% Compression
    % Don't expand matrix if result matrix is smaller
    %nx = min(nx, N(i));
    %nt = min(nt, M);
%     [Z, xcomp, tcomp] = compress(x, t, U, nx, nt);
%     result(:,:,i) = Z;
end
%% Plotting

% legenddata = cellstr(num2str((N)'));
% ylabel('Relative fault')
% xlabel('Time')
% legend(legenddata,'Location','NorthWest')

% figure
% hold on
% for i = 1:nM
%     plot(xcomp,result(end,:,i)','color',color(i,:))
% end
% plot(xcomp,ref(xcomp,T),'k')
% legenddata(nM + 1 ) = cellstr('Analytical solution');
% legend(legenddata,'Location','NorthWest')
% xlabel('Travelling peakon')

p = polyfit(log(h), log(errors), 1)

figure
hold on
loglog(h, errors,'b*-')

%hold on
% plot(h,fliplr(h) - max(h) + max(errors),'r*-')
% plot(h,fliplr(h.^2) - max(h.^2) + max(errors),'g*-')
% legend('Convergence of our method', 'Linear convergence','Quadratic convergence','Location','NorthEast')

