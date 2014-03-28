function [ ] = animatedplot(x, t, U)
% Target framerate
FRAMERATE = 30;
FRAMELENGTH = 1 / FRAMERATE;

% Plot initial data
figure
plothandle = plot(x, U(1, :));
Umax = max(U(1, :));
Umin = min(min(U));
ylim([Umin, Umax]);
for j = 2:length(t)
    set(plothandle, 'YData', U(j, :))
    drawnow
    pause(FRAMELENGTH)
end
end