clear all; close all
%% The purpose of this script is to calculate channel length in the control and treatment experiments. 
% These data are described in Table 1 and plotted in Figure 3b in Sanks et al.
% (2023) submitted to Earth Surface Dynamics
% We will caclculate the channel length by determining the distance from
% the apex to the furthest channel pixel

%% First, lets set our relative path
cd '../' %We should have opened this code from 'C:\Users\XXXXXXX\Documents\GitHub\TDWB_19_2_Channels\code'

%% Load data
cd '.\data'
load('CM_18.mat');
load('CM_19.mat');
cd '..\code'

%% Define and set parameters
% control
nt_18 = size(CM_18,3); % number of time steps in data set
xentrance_18 = 109; % x grid node location of the entrance channel
yentrance_18 = 271; % y grid node location of the entrance channel
CM_18(CM_18 == 0) = NaN; % all data outside channels is now NaN

% treatment
nt_19 = size(CM_19,3); % number of time steps in data set
xentrance_19 = 214; % x grid node location of the entrance channel (x is down dip)
yentrance_19 = 397; % y grid node location of the entrance channel (y is strike)
CM_19(CM_19 == 0) = NaN; % all data outside channels is now NaN

%% Calculate channel length by finding the distance from the apex to the furthest channel pixel
% control
Lc_18 = []; % initialize matrix
for i = 1:nt_18 % loop through time
    chan = CM_18(:,:,i); % channel map at time i
    dist_18 = [];
    for xx = 1:length(chan(:,1)) % loop through x distances 
        for yy = 1:length(chan(1,:)) % loop through y distances
            if ~isnan(chan(xx,yy)) % if channel pixel
            dist = 5*(sqrt((abs(xentrance_18-xx))^2 + abs((yentrance_18-yy))^2)); % calculate distance from apex
            dist_18 = [dist_18;dist]; % save distances
            end
        end
    end
    maxdist = (max(dist_18))/10^3; % find max distance from apex for each timestep (meters)
    Lc_18 = [Lc_18;maxdist]; % channel length is approximated as max distance form apex
end 

% treatment
Lc_19 = []; % initialize matrix
for i = 1:nt_19 % loop through time
    chan = CM_19(:,:,i); % channel map at time i
    dist_19 = [];
    for xx = 1:length(chan(:,1)) % loop through x distances    
        for yy = 1:length(chan(1,:)) % loop through y distances
            if ~isnan(chan(xx,yy)) % if channel pixe;
                dist = 5*(sqrt((abs(xentrance_19-xx))^2 + abs((yentrance_19-yy))^2)); % calculate distance from apex
                dist_19 = [dist_19;dist]; % save distances
            end
        end
    end
    if isempty(dist_19) % if no channel map
        maxdist = NaN; % no channel length
        Lc_19 = [Lc_19;maxdist];
    else
        maxdist = (max(dist_19))/10^3; % find max distance from apex for each timestep (meters)
        Lc_19 = [Lc_19;maxdist]; % channel length is approximated as max distance form apex
    end
end 

%% Calculate statistics (mean and standard deviation of Lc) for Table 1
mean_Lc_18 = mean(Lc_18, 'omitnan');
mean_Lc_19 = mean(Lc_19, 'omitnan');
stdev_Lc_18 = std(Lc_18, 'omitnan');
stdev_Lc_19 = std(Lc_19, 'omitnan');

%% Plot the data
% data for violin plot of channel length (Lc)
G = [ones(size(Lc_18)); 2*ones(size(Lc_19))];
X = [Lc_18; Lc_19];

% Violin plot for channel length experiments
% this comes from 10.5281/zenodo.4559847
fig = figure();
violinplot(X,G);
ylabel('channel length (m)');
ylim([0 3])
set(gcf, 'PaperUnits', 'inches');
y_width=7.25 ;x_width=9.125
set(gcf, 'PaperPosition', [0 0 x_width y_width]);
saveas(gcf,'../figures/esurf_Figure3b.pdf')
