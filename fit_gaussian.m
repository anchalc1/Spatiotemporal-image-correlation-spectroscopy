function [fitParams,flag] = fit_gaussian(gData,coords)
%  fit_gaussian finds the best fit 2D Gaussian fit.
%
% fit_gaussian uses lsqcurvefit to find a Gaussian fit for GDATA plotted
% versus coordinates COORDS. 
% GDATA is a 2D matrix
% COORDS = [X,Y] X and Y are 2D arrays of the type returned by meshgrid
%
% FITPARAMS = = [a,sigma,b,xd,yd]
% FLAG is the flag returned by lsqcurvefit. See lsqcurvefit documentation.
%
% Copyright (c) Molly J. Rossow 2014

X = coords(:,:,1);
Y = coords(:,:,2);


% Guess that the gaussian peak is located at the maximum
[~,ind] = max(gData(:));
xd = X(ind);
yd = Y(ind);
initParams = [max(gData(:)),5,xd,yd,min(gData(:))];
lb = [0,0,0,min(X(:)),min(Y(:))];
ub = [10,10,inf,max(X(:)),max(Y(:))];
options = optimoptions('lsqcurvefit', 'TolFun',1e-16, ...
    'MaxFunEvals', 10000);
[fitParams, ~, ~, flag] = lsqcurvefit(@gaussian_2D,initParams,...
    coords,gData,lb,ub,options);
