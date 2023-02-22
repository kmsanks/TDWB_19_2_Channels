clear all; close all

%% The purpose of this script is to calculate the fraction of total, overbank, and channelized flow on the delta top
% These data are plotted in Figure 1 and Table 1 in Sanks et al.
% (2023) submitted to Earth Surface Dynamics
% total flow is all the "wet" pixels above sea level
% channelized flow is the flow contained within the channel maps
% overbank flow is total flow - channelized flow

%% First, lets set our relative path
cd '../' %We should have opened this code from 'C:\Users\XXXXXXX\Documents\GitHub\TDWB_19_2_Channels\code'

%% Load data
% control
cd './data'
load('ZD_18.mat'); %topography; elevation data (mm)
load('CM_18.mat'); %channel maps
load('flowscreen18.mat'); %area covered by flow

% treatment
load('ZD_19.mat'); %this will be used to create a binary for the basin
load('ZW_19.mat'); %topography; elevation data (mm) from the wet scans 
load('CM_19.mat'); %channel maps
load('flowscreen19.mat'); %area covered by flow
cd '../code'

%% Set parameters
% control 
nx_18 = size(ZD_18, 1); %number of x locations on map
ny_18 = size(ZD_18,2); %number of y locations on map
nt_18 = size(ZD_18,3); %number of time steps in data set
dt_18 = 1; %delta t of time steps (hr)

% treatment
nx_19 = size(ZW_19,1); %number of x locations on map
ny_19 = size(ZW_19,2); %number of y locations on map
nt_19 = size(ZW_19,3); %number of time steps in data set
dt_19 = 1; %delta t of time steps for the wet scans (hr)

% both
dx = 5; %5 mm grid cells in x
dy = 5; %5 mm grid cells in y
baselevel_rr = 0.25; %base level rise rate (mm/hr)
ocean_zero = 25; %ocean elevation at beginning of experiment (mm)

%% Elevation mask, so flow can be analyzed only on the area about sea level
% We will use a boundary method, so floating mats are excluded from
% subsequent analyses

% control
z18 = []; % initialize elevation matrix
terr_area18 = []; % initialize terrestrial area binary matrix
for i= 1:nt_18; % loop through all timesteps
    %What is sea level at time i
    sl = (i*0.25*dt_18)+25; %first wet scan starts 48 minutes into the hour so i = i is a good approximation
    elevationmask = ZD_18(:,:,i); 
    elevationmask(elevationmask == 0) = NaN;
    elevationmask_rslr = elevationmask - sl;
    elevationmask_rslr(elevationmask_rslr < 0) = NaN;
    z18(:,:,i) = elevationmask_rslr; %elevation relative to sea level

    % Now make the shoreline boundary that removes any floating mats from
    b = elevationmask_rslr; 
    
    % Use this for shoreline boundary
    b(b > 0) = 1;
    b(b <= 0) = 0;  
    
    b(isnan(b)) = 0;% address the nan values that become a problem later when creating a binary mask
    B = bwboundaries(b,'noholes'); % Find the boundary of the matrix...creates n cell arrays of varying sizes
    [msize, mindex] = max(cellfun('size',B,1)); % Find the cell array with the largest size...this is the cell array that contains the shoreline locations
    C = B{mindex,1}; % Get the length of the largest cell array that contains the shoreline locations

    C2 = fliplr(C);
    imagesc(b(:,:,1))
    terr = drawpolygon('Position',C2);
    terr_area18(:,:,i) = createMask(terr);
end

% treatment
z19 = []; % initialize elevation matrix
terr_area19 = []; % initialize terrestrial area binary matrix
for i= 1:nt_19; % loop through all timesteps
    %What is sea level at time i
    sl = (i*0.25*dt_19)+25; %first wet scan starts 48 minutes into the hour so i = i is a good approximation
    elevationmask = ZW_19(:,:,i); 
    elevationmask(elevationmask == 0) = NaN;
    elevationmask_rslr = elevationmask - sl;
    elevationmask_rslr(elevationmask_rslr < 0) = NaN;
    z19(:,:,i) = elevationmask_rslr; %elevation relative to sea level

    % Now make the shoreline boundary that removes any floating mats from
    b = elevationmask_rslr; 
    
    % Use this for shoreline boundary
    b(b > 0) = 1;
    b(b <= 0) = 0;  
    
    b(isnan(b)) = 0;% address the nan values that become a problem later when creating a binary mask
    B = bwboundaries(b,'noholes'); % Find the boundary of the matrix...creates n cell arrays of varying sizes
    [msize, mindex] = max(cellfun('size',B,1)); % Find the cell array with the largest size...this is the cell array that contains the shoreline locations
    C = B{mindex,1}; % Get the length of the largest cell array that contains the shoreline locations

    C2 = fliplr(C);
    imagesc(b(:,:,1))
    terr = drawpolygon('Position',C2);
    terr_area19(:,:,i) = createMask(terr);
end

%% Calculate the flow (total, channel, and overbank) properties 
% 1 is flow and 0 is no flow

% control 
flowscreen18 = flowscreen18.*z18.*terr_area18; % multiply flow, elevation relative to sea level, and terrestrial delta
flowscreen18(flowscreen18 > 0) = 1; %1 is flow, 0 is no flow, NaN is not on terrestrial delta

% initialize empty matricies
chanarea18 = [];
deltaarea18 = [];
flowarea18 = [];
obarea18 = [];
for i =1:nt_18
  chan = CM_18(:,:,i);
  flow = flowscreen18(:,:,i);
  flowtot = flow + chan; %deep parts of channel are sometimes not included in flow
  flowtot(flowtot>=1) = 1;
  ob = flowtot-chan; %0 channel or no flow, 1 ob flow, NaN outside basin or not terrestrial
  area = terr_area18(:,:,i)+chan; % add channels since backwater length not captures in terrestrial area
  area(area>=1)=1; 
  chanarea18(i) = sum(chan(:), 'omitnan')*(2.5*10^-5); %m^2
  obarea18(i) = sum(ob(:), 'omitnan')*(2.5*10^-5); %m^2
  flowarea18(i) = sum(flowtot(:), 'omitnan')*(2.5*10^-5); %m^2
  deltaarea18(i) = sum(area(:), 'omitnan')*(2.5*10^-5); %m^2
end

% treatment 
% 1 is flow (and outside basin) and 0 is no flow
flowscreen19 = flowscreen19.*z19.*terr_area19; % multiply flow, elevation relative to sea level, and terrestrial delta
flowscreen19(flowscreen19 > 0) = 1; %1 is flow, 0 is no flow, NaN is not on terrestrial delta

% initialize empty matrices
chanarea19 = [];
deltaarea19 = [];
flowarea19 = [];
obarea19 = [];
%Visual spot checking shows we only need to remove the following 76 hours due to issues with cart messing up flow area:
idx = [79,85,87,89,91,93,95,97,99,...
    101,105,107,109,111,113,121,123,125,127,131,133,135,139,143,147,151,153,155,158,161,163,165,167,171,175,183,187,191,192,195,199,...
    203,207,215,219,227,230,231,239,243,247,248,251,255,259,267,271,291,293,294,297,...
    301,309,310,349,381,385,389,390,...
    401,405,417,443,449,450,...
    541];
for i =1:nt_19
  chan = CM_19(:,:,i);
  area = terr_area19(:,:,i) + chan; % add channels because backwater not included in terrestrial map
  area(area>=1)=1;
  % skip time steps with no channel map or cart impaired flow
  if sum(chan(:), 'omitnan') == 0 | ismember(i,idx) %if missing a channel map or cart impaired flow map, flow stats should be NaN
     chanarea19(:,i) = NaN;
     deltaarea19(:,i) = length(area(~isnan(area)))*(2.5*10^-5); %m2
     flowarea19(:,i) = NaN;
     obarea19(:,i) = NaN;
  else % calculate metrics for all other timesteps
      flow = flowscreen19(:,:,i);
      flowtot = flow + chan; %deep parts of channel are sometimes not included in flow
      flowtot(flowtot>=1) = 1;
      ob = flowtot-chan; %0 channel or no flow, 1 ob flow, NaN outside basin or not terrestrial
      %we need to ignore nans from the cart here CHECK
      %tmp = cat(3,flow,chan);  CHECK
      %flowmap19_tot = sum(tmp,3, 'omitnan'); CHECK
      %flowmap19_tot = flow + chan; %deep parts of channel are sometimes not included in flow
      %flowmap19_tot(flowmap19_tot>=1) = 1; CHECK
      chanarea19(i) = sum(chan(:), 'omitnan')*(2.5*10^-5); %m^2
      obarea19(i) = sum(ob(:), 'omitnan')*(2.5*10^-5); %m^2
      flowarea19(i) = sum(flowtot(:), 'omitnan')*(2.5*10^-5); %m^2
      deltaarea19(i) = sum(area(:), 'omitnan')*(2.5*10^-5); %m^2
  end
end


%% Calculate fraction of flow

% control
flowfrac18 = flowarea18./deltaarea18; %fraction of total flow on delta top
chanfrac18 = chanarea18./deltaarea18; %fraction of channel flow delta top
obfrac18 = obarea18./deltaarea18; %fraction of overbank flow on delta top
chan_overbank_frac18 = chanarea18./obarea18; %ratio of channel flow to overbank flow

%treatment
flowfrac19 = flowarea19./deltaarea19; %fraction of total flow on delta top
chanfrac19 = chanarea19./deltaarea19; %fraction of channel flow delta top
obfrac19 = obarea19./deltaarea19; %fraction of overbank flow on delta top
chan_overbank_frac19 = chanarea19./obarea19; %ratio of channel flow to overbank flow

%% Calculate mean statistics
% area
% channel area and standard deviation
mean_chan_area18 = mean(chanarea18, 'omitnan');
stdev_chan_area_18 = std(chanarea18, 'omitnan');
mean_chan_area19 = mean(chanarea19, 'omitnan');
stdev_chan_area_19 = std(chanarea19, 'omitnan');
% overbank area and standard deviation
mean_ob_area18 = mean(obarea18, 'omitnan');
stdev_ob_area18 = std(obarea18, 'omitnan');
mean_ob_area19 = mean(obarea19, 'omitnan');
stdev_ob_area19 = std(obarea19, 'omitnan');
% total flow area and standard deviation
mean_tot_area18 = mean(flowarea18, 'omitnan');
stdev_tot_area18 = std(flowarea18, 'omitnan');
mean_tot_area19 = mean(flowarea19, 'omitnan');
stdev_tot_area19 = std(flowarea19, 'omitnan');

% fraction
% channel frac and standard deviation
mean_chan_frac18 = mean(chanfrac18, 'omitnan');
stdev_chan_frac_18 = std(chanfrac18, 'omitnan');
mean_chan_frac19 = mean(chanfrac19, 'omitnan');
stdev_chan_frac_19 = std(chanfrac19, 'omitnan');
% overbank frac and standard deviation
mean_ob_frac18 = mean(obfrac18, 'omitnan');
stdev_ob_frac18 = std(obfrac18, 'omitnan');
mean_ob_frac19 = mean(obfrac19, 'omitnan');
stdev_ob_frac19 = std(obfrac19, 'omitnan');
% total flow frac and standard deviation
mean_tot_frac18 = mean(flowfrac18, 'omitnan');
stdev_tot_frac18 = std(flowfrac18, 'omitnan');
mean_tot_frac19 = mean(flowfrac19, 'omitnan');
stdev_tot_frac19 = std(flowfrac19, 'omitnan');

% ratio of channel to overbank flow
mean_ratio18 = mean(chan_overbank_frac18, 'omitnan');
std_ratio18 = std(chan_overbank_frac18, 'omitnan');
mean_ratio19 = mean(chan_overbank_frac19, 'omitnan');
std_ratio19 = std(chan_overbank_frac19, 'omitnan');

%% Find a control and treatment timestep that are representative to average to show above bottom plot
idx18 = find(abs(chanfrac18-mean_chan_frac18)<=0.0005 & min(abs(obfrac18-mean_ob_frac18)) <= 0.0005);
idx19 = find(abs(chanfrac19-mean_chan_frac19)<=0.0005 & min(abs(obfrac19-mean_ob_frac19)) <= 0.0005);

% results:
% control = 83,102,181,414,451,453,468,535
% treatment = 68, 360, 400

%% Plot the data
% data for violin plot of channel and overbank fraction
G = [ones(size(chanfrac18)), 2*ones(size(chanfrac19)), 3*ones(size(obfrac18)), 4*ones(size(obfrac19))];
X = [chanfrac18, chanfrac19, obfrac18, obfrac19];

% Violin plot for channel and overbank flow fraction for both experiments
% this comes from 10.5281/zenodo.4559847

fig = figure();
violinplot(X,G);
ylabel('fraction of delta top covered with flow (-)');
set(gcf, 'PaperUnits', 'inches');
x_width=7.25; 
y_width=9.125;
set(gcf, 'PaperPosition', [0 0 x_width y_width]);
saveas(gcf,'../figures/esurf_Figure1.pdf') % this was modified in illustrator for visual purposes

%% Statistical tests
% is the difference in flow fraction statistically significant?
[h, p, ci, stats] = ttest2(flowfrac18, flowfrac19)
[h_un, p_un, ci_un, stats_un] = ttest2(flowfrac18, flowfrac19, 'Vartype', 'unequal')
