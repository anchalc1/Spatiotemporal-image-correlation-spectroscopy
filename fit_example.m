% Example of fitting a single frame of STICS data using lsqcurvefit to find
% the offset of the STICS peak. The offset corresponds to the distance
% traveled between frames.
%
%  Copyright (c) Molly Rossow 2014

% Create fake data
s = 128; % Image size
halfS = round(s/2);
x = -halfS:halfS;
y = -halfS:halfS;
[X,Y] = meshgrid(x,y); % X and Y for input to gaussian_2D
coords = cat(3,X,Y);
params = [10,5,20,15,25]; %[a,sigma,b,xd,yd]
% gRaw with random noise
gData = gaussian_2D(params,coords) + 10*rand(size(X));


% Find best fit parameters
[fitParams,flag] = fit_gaussian(gData,coords);

% Calculate 2D gaussian with fit parameters. 
gFit = gaussian_2D(fitParams,coords);

% Dispaly raw data
figure
subplot(2,1,1)
imagesc(gData)
title('G Data')

% Dispaly fit 
subplot(2,1,2)
imagesc(gFit)
title('G Fit')
