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
I = 8:2:14;
h = zeros(1, length(I));
N = 2.^I;
ts = 2.^(-I);
% Maximum time value
T = 0.01;
xmin = -10;
xmax = 35;

% Compression settings
% nx = number of x values in compressed matrix
nx = 800;
% nt = number of t values in compressed matrix
nt = 1;

%% Initial condition and reference solution
ref = @(x,t) exp(-abs(x-t));
initial = @(x) ref(x, 0);

figure
hold on
color = hsv(length(I));

%% Loop to solve equation and check convergence
errors = zeros(1, length(I));
finalvalue = zeros(length(I), nx);
normsplot = zeros(length(I), nx);
for j = 1:length(I);
    [U, x, t] = holdenraynaud(N(j), T, [xmin, xmax], initial);
    
    h(j) = x(2) - x(1);
    %generate reference solution:
    norms = zeros(1,length(t));
    for l = 1:length(t)
       norms(l) = sqrt(h(j))*norm(ref(x, t(l)) - U(l,:),2);
    end
    plot(t,norms,'color',color(j,:))
    
    % Error function norm
    %errors(j) = sqrt(h(j))*norm(ref(x,T)-U(end,:));
    errors(j) = norms(end);
    
    %% Compression of the function values at time T
    [Z, xcomp, tcomp] = compress(x, t, U, nx, nt);
    finalvalue(j, :) = Z;   
    
end

p = polyfit(log(h), log(errors), 1)

%% Plotting
% info about normplots over time
legenddata = cellstr(num2str((N)'));
ylabel('Relative fault')
xlabel('Time')
legend(legenddata,'Location','NorthWest')
%print(gcf,'-dpng','-r400','erroroftime')

% plot of convergence
figure
plot(log(h),log(errors),'b*')
hold on
x = linspace(log(h(1)),log(h(end)));
f = p(1)*x + p(2);
f2 = 1*x + p(2);
plot(x,f,'b')
plot(x,f2,'r')
xlabel('log of stepsize h')
ylabel('log of error')
%print(gcf,'-dpng','-r400','loglog')

% plot peakons at time T
legenddata(length(I) + 1 ) = cellstr('Analytical solution');
figure
hold on
for i = 1:length(I)
    plot(xcomp,finalvalue(i,:),'color', color(i,:),'LineWidth',1.3)
end
% plot analytical solution
plot(xcomp,ref(xcomp,T),'k','LineWidth',1.3)
legend(legenddata,'location', 'best')
%print(gcf,'-dpng','-r400','attimeT')

