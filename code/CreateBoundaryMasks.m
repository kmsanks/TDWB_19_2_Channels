% This script will calculate compensation index
% This script will plot transects that show how the shoreline changes and
% will use FFT to look at any signals within these time series 
% Current version 1/31/22
clear all
close all
%% Import data
% mac

load('/Users/josesilvestre/OneDrive - Tulane University/numerical techniques in strat/ZD_19_2_dry.mat')
ZD = ZD_19_2_dry; 

%% Subtract elevation to find true ocean level

for i = 1:size(ZD,3)
    test = ZD(:,:,i); 
    test(test == 0) = nan; 
    newZD(:,:,i) = test-(0.25*(i*2)); 
end

%% Convert topo data to delta area mask


% Determine terrestrial and ocean locations for every time step
% 1's = terrestrial
% 0's = Ocean
% Note: These 1's and 0's are not binary since we will need to add the
% value at each pixel through time

for i=1:size(ZD_18,3)
    i
    newB = [];
    test2 = z18(:,:,i); 
    
%   Use this for shoreline boundary
    test2(test2 > 0) = 1;
    test2(test2 <= 0) = 0;  
    
    % Use this for upper marsh boundary
%     test2(test2 <= 5) = 0;
%     test2(test2 > 5) = 1;  
%     
%     % Use this for lower marsh boundary
%     test2(test2 > -9) = 1;
%     test2(test2 <= -9) = 0;  

    test2(isnan(test2)) = 0;            % address the nan values that become a problem later when creating a binary mask
    newZD2(:,:,i) = test2;
    B = bwboundaries(test2,'noholes'); % Find the boundary of the matrix...creates n cell arrays of varying sizes
    [msize, mindex] = max(cellfun('size',B,1)); % Find the cell array with the largest size...this is the cell array that contains the shoreline locations
    C = B{mindex,1}; % Get the length of the largest cell array that contains the shoreline locations

    C2 = fliplr(C);
    imagesc(test2(:,:,1))
    h = drawpolygon('Position',C2);
    h2(:,:,i) = createMask(h);
end
 
