function [ V, x, y ] = compress( X, Y, U, xnum, ynum )
%% COMPRESS(U, xnum ynum) 
%       Compresses the numerical grid matrix U into a matrix V of size 
%       xnum x ynum elements. Assumes X and Y are linearly spaced.

x = linspace(min(X), max(X), xnum);
y = linspace(min(Y), max(Y), ynum);

[xq, yq] = meshgrid(x, y);
V = interp2(X, Y, U, xq, yq);

assert(size(V, 1) == ynum);
assert(size(V, 2) == xnum);

end