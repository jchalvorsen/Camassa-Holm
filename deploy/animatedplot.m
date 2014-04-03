function [ ] = animatedplot(x, t, varargin)
% Target framerate
FRAMERATE = 30;
FRAMELENGTH = 1 / FRAMERATE;

% Every cell in varargin is a matrix of (x, t) data
% Number of data sets
N = length(varargin);

% Number of time points
T = 0;
% Let U0 hold the data initial data across all series
U0 = zeros(N, length(x));

% Max/min values for axis bounds
Umax = -Inf;
Umin = Inf;

for j = 1:N
    T = max(T, size(varargin{j}, 1));
    U0(j, :) = varargin{j}(1, :);
    Umax = max(Umax, max(max(varargin{j})));
    Umin = min(Umin, min(min(varargin{j})));
end

% Plot initial data and proceed to animate the rest
figure
plothandle = plot(x, U0);
ylim([Umin, Umax]);
for j = 2:T
    for k = 1:N
        set(plothandle(k), 'YData', varargin{k}(j, :));
    end
    
    drawnow
    pause(FRAMELENGTH)
end

end