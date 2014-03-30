function [ ] = animatedplot(x, t, varargin)
% Target framerate
FRAMERATE = 30;
FRAMELENGTH = 1 / FRAMERATE;

% Every cell in varargin is a matrix of (x, t) data
% Number of data sets
N = length(varargin);

% Number of time points
T = size(varargin{1}, 1);
for j = 2:N
    T = max(T, size(varargin{j}, 1));
end

% Let U hold the data for each time step across all series
U = zeros(N, length(x));
for j = 1:N
    U(j, :) = varargin{j}(1, :);
end

% Plot initial data
figure
plothandle = plot(x, U);
%Umax = max(U(1, :));
Umax = max(max(U));
Umin = min(min(U));
ylim([Umin, Umax]);
for j = 2:T
    for k = 1:N
        set(plothandle(k), 'YData', varargin{k}(j, :));
    end
    
    drawnow
    pause(FRAMELENGTH)
end

end