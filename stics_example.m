% Example STICS analysis using simulated data. 
%
% copyright Molly J. Rossow 2014

% Load simulated data
load('simulatedData.mat','data');

% Select a subregion of the data
subRegion = data(100:201,100:201,:);

% STICS calucation
sticsG = calc_stics(subRegion,6,1);

% Plot results
figure
for t = 1:6
   subplot(3,2,t)
   imagesc(sticsG(:,:,t))
   title(['shift = ', num2str(t)])
   xlabel('\xi')
   ylabel('\eta')
end
