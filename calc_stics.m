function sticsG = calc_stics(matrix3Dim,shifts,flag)
% calc_stics calculates the spatial temporal image correlation spectroscopy
% for a series of frames. 
%
% STICSG - calc_stics(MATRIX3DIM, SHIFTS, FLAG) 
% MATRIX3DIM - 3 dimensional matrix, [y,x,time] (x and y must be even and
% equal to each other)
% SHIFTS - integer, number of STICS shifts to calculate 
% FLAG  = 0, no interpolation i.e. sticsG is [x/2,y/2,shifts]. flag = 1
% sticsG is interpolated at the center line for symmetrical plotting.
% 
% Background subtraction is performed before STICS analysis. 
%
% copyright Molly J. Rossow 2007 

if nargin < 3
    flag = 0;
end

s = size(matrix3Dim);
xdim = s(2);
ydim = s(1);
numFrames = s(3);

% Check image size
if ~(xdim == ydim)
end

%Normalizes all the frames by substracting the average of all frames.
%Prevents frames from having negative values by keeping adding in the
%average of each frame. 
averageAllFrame = mean(matrix3Dim,3);
averageEachFrame = mean(mean(matrix3Dim,1),2);
matrix3Dim = double(matrix3Dim);
matrix3Dim = matrix3Dim - repmat(averageAllFrame,[1,1,numFrames]) + ...
    repmat(averageEachFrame,[ydim,xdim,1]);

fTransform = fft(matrix3Dim); %one dimensionaly fft
sumOfMatrix = sum(sum(matrix3Dim,1),2);
clear matrix3Dim;
repsums = repmat(sumOfMatrix,[ydim,xdim]); 
clear sumOfMatrix

sticsG = zeros(ydim, xdim, shifts);

for shift = 1:shifts
    P = permute(fft(permute(fTransform(:,:,shift:end),[2,1,3])),[2,1,3])...
        .*conj(permute(fft(permute(fTransform(:,:,1:end-shift+1),...
        [2,1,3])),[2,1,3]));
    gRaw = permute(ifft(permute(ifft(P,'nonsymmetric'),[2,1,3]),...
        'nonsymmetric'),[2,1,3]);
    clear P
    gRaw = gRaw(1:ydim,1:xdim,:);
    gNorm = gRaw.*xdim^2./(repsums(:,:,shift:end).*...
        repsums(:,:,1:end-shift+1))-1;
    temp = mean(gNorm,3);
    sticsG(:,:,shift) = [temp(ydim/2+1:end,xdim/2+1:end), ...
        temp(ydim/2+1:end,1:xdim/2); temp(1:ydim/2,xdim/2+1:end), ...
        temp(1:ydim/2,1:xdim/2)];
end
sticsG = real(sticsG); % Eliminates rounding error. The imaginary part 
                       % should be very small.

if flag == 1
    % Interpoloate 0 shift points
    s = size(sticsG);
    newSticsG = zeros(s(1)+1,s(2)+1,s(3));
    newSticsG(1:s(1)/2,1:s(2)/2,:) = ...
        sticsG(1:s(1)/2,1:s(2)/2,:); %upper left quadrant
    newSticsG(1:s(1)/2,s(2)/2+2:end,:) = ...
        sticsG(1:s(1)/2,s(2)/2+1:end,:); %upper right quadrant
    newSticsG(s(1)/2+2:end,1:s(2)/2,:) = ...
        sticsG(s(1)/2+1:end,1:s(2)/2,:); %lower left quadrant
    newSticsG(s(1)/2+2:end,s(2)/2+2:end,:) = ...
        sticsG(s(1)/2+1:end,s(2)/2+1:end,:); %lower left quadrant
    %top vertical line
    newSticsG(1:s(1)/2,s(2)/2+1,:) = ...
        (sticsG(1:s(1)/2,s(2)/2,:)+sticsG(1:s(2)/2,s(2)/2+1,:))./2; 
    %bottom vertical line
    newSticsG(s(1)/2+2:end,s(2)/2+1,:) = ...
        (sticsG(s(1)/2+1:end,s(2)/2,:)+sticsG(s(2)/2+1:end,s(2)/2+1,:))./2; 
    %right horizotal
    newSticsG(s(1)/2+1,s(2)/2+2:end,:) = ...
        (sticsG(s(1)/2,s(2)/2+1:end,:)+sticsG(s(2)/2,s(2)/2+1:end,:))./2;
     %left horizontal
    newSticsG(s(1)/2+1,1:s(2)/2,:) = ...
        (sticsG(s(1)/2,1:s(2)/2,:)+sticsG(s(2)/2,1:s(2)/2,:))./2;

    newSticsG(s(1)/2+1,s(2)/2+1,:) = (sticsG(s(1)/2,s(2)/2,:) + ...
        sticsG(s(1)/2+1,s(2)/2,:) + sticsG(s(1)/2,s(2)/2+1,:) + ...
        sticsG(s(1)/2+1,s(2)/2+1,:))/4;
    sticsG = newSticsG;
end
