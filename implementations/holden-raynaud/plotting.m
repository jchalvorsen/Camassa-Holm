close all
clear all

% plotting analytical solution of the equation

peakon = @(x,t,c,h) c*exp(-abs(x - h - c*t));
antipeakon = @(x,t,c,h) c*exp(-abs(x - h + c*t));
doublePeakon = @(x,t) peakon(x,t,1,-1) + peakon(x,t,1,1);
peakonOvertake = @(x,t) peakon(x,t,2,0) + peakon(x,t,1,5);
peakonAntipeakon = @(x,t) peakon(x,t,1,0) - antipeakon(x,t,1,20);

xmin = 0;
xmax = 10;
x = linspace(xmin,xmax,2000);

% Single peakon
t = [ 0 5 10 ];
p1 = figure;
hold on
plot(x,peakon(x,t(1),1,0),'b')
plot(x,peakon(x,t(2),1,0),'b--')
%legend('t = 0', 't = 5', 'Location', 'NorthEast')


% Double peakon
p2 = figure;
hold on
plot(x,doublePeakon(x,t(1)),'b')
plot(x,doublePeakon(x,t(2)),'b--')
%legend('t = 0', 't = 5', 'Location', 'NorthEast')


% Double peakon interaction
xmin = 0;
xmax = 25;
x = linspace(xmin,xmax,2000);
p3 = figure;
hold on
plot(x,peakonOvertake(x,t(1)),'b')
plot(x,peakonOvertake(x,t(2)),'b--')
plot(x,peakonOvertake(x,t(3)),'b-.')
%legend('t = 0', 't = 5', 't = 10', 'Location', 'NorthEast')


t = [0 10 15];
% Peakon antipeakon interaction.
p4 = figure;
subplot(3,1,1)
plot(x,peakonAntipeakon(x,t(1)),'b')
%legend('t = 0', 'Location', 'NorthEast')
subplot(3,1,2)
plot(x,peakonAntipeakon(x,t(2)),'b--')
axis([xmin xmax -1 1])
%legend('t = 10', 'Location', 'NorthEast')
subplot(3,1,3)
plot(x,peakonAntipeakon(x,t(3)),'b-.')
%legend('t = 15', 'Location', 'NorthEast')


% print images to file
% print(p1,'-dpng','-r400','peakon')
% print(p2,'-dpng','-r400','doublepeakon')
% print(p3,'-dpng','-r400','peakonovertake')
% print(p4,'-dpng','-r400','peakonantipeakon')


% plot area under curve for peakon and another initial condition
p5 = figure;
hold on
exps = exp(-((x-(xmax-xmin)/2)/4).^2);
area(x,exps)
areaundersine = 7.09;
c = areaundersine / 2;
axis([xmin xmax 0 c])
p6 = figure;
area(x,peakon(x,0,c,(xmax-xmin)/2))
axis([xmin xmax 0 c])

print(p5,'-dpng','-r400','areainitial')
print(p6,'-dpng','-r400','areapeakon')


