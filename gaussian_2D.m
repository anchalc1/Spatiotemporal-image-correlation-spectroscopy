function gfit = gaussian_2D(params, coords)
% gaussian_2D 2D Gaussian function calculated using the parameters PARAMS
% and the coordinates COORDS
% 
% params = [a,sigma,b,xd,yd]
% coords = [X,Y] X and Y are 2D arrays of the type returned by meshgrid
%
% Copyright (c) Molly J. Rossow 2014

a = params(1);
sigma = params(2);
b = params(3);
xd = params(4); % x offset
yd = params(5); % y offset

X = coords(:,:,1);
Y = coords(:,:,2);

gfit = a * exp(-1 * (((X-xd).^2/sigma^2) + ((Y-yd).^2/sigma^2))) + b;


end