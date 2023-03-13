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
load('area18.mat')
% treatment
load('ZD_19.mat'); %this will be used to create a binary for the basin
load('ZW_19.mat'); %topography; elevation data (mm) from the wet scans 
load('CM_19.mat'); %channel maps
load('flowscreen19.mat'); %area covered by flow
load('area19.mat')
cd '../code'

%% Set parameters
% control 
nx_18 = size(ZD_18, 1); %number of x locations on map
ny_18 = size(ZD_18,2); %number of y locations on map
nt_18 = size(ZD_18,3); %number of time steps in data set
dt_18 = 1; %delta t of time steps (hr)
xentrance_18 = 109; % x grid node location of the apex
yentrance_18 = 271; % y grid node location of the apex

% treatment
nx_19 = size(ZW_19,1); %number of x locations on map
ny_19 = size(ZW_19,2); %number of y locations on map
nt_19 = size(ZW_19,3); %number of time steps in data set
dt_19 = 1; %delta t of time steps for the wet scans (hr)
xentrance_19 = 214; % x grid node location of the apex (x is down dip)
yentrance_19 = 397; % y grid node location of the apex (y is strike)

% both
dx = 5; %5 mm grid cells in x
dy = 5; %5 mm grid cells in y
baselevel_rr = 0.25; %base level rise rate (mm/hr)
ocean_zero = 25; %ocean elevation at beginning of experiment (mm)

%% get boundary of basin, so we don't have issues later
% control
basin18 = ZD_18(:,:,1);
basin18(basin18 == 0) = NaN; % everything outside basin is NaN
basin18(~isnan(basin18)) = 1;% everything inside basin is 1

% treatment
basin19 = ZD_19(:,:,1);
basin19(basin19 == 0) = NaN; % well shadow is NaN
basin19(~isnan(basin19)) = 1;% everything inside basin is 1 except for well

%% Elevation mask, so flow can be analyzed only on the area about sea level
% We will use a boundary method, so floating mats are excluded from
% subsequent analyses

% % control
% area18 = []; % initialize terrestrial area binary matrix
% for i= 1:nt_18 % loop through all timesteps
%     %What is sea level at time i
%     i
%     sl = (i*0.25*dt_18)+25; %first wet scan starts 48 minutes into the hour so i = i is a good approximation
%     elevationmask = ZD_18(:,:,i); 
%     elevationmask(elevationmask == 0) = NaN;
%     elevationmask_rslr = elevationmask - sl;
%     elevationmask_rslr(elevationmask_rslr < 0) = NaN;
% 
%     % Now make the shoreline boundary that removes any floating mats from
%     b = elevationmask_rslr; 
%     
%     % Use this for shoreline boundary
%     b(b > 0) = 1;
%     b(b <= 0) = 0;  
%     
%     b(isnan(b)) = 0;% address the nan values that become a problem later when creating a binary mask
%     B = bwboundaries(b,'noholes'); % Find the boundary of the matrix...creates n cell arrays of varying sizes
%     [msize, mindex] = max(cellfun('size',B,1)); % Find the cell array with the largest size...this is the cell array that contains the shoreline locations
%     C = B{mindex,1}; % Get the length of the largest cell array that contains the shoreline locations
% 
%     C2 = fliplr(C);
%     imagesc(b(:,:,1))
%     terr = drawpolygon('Position',C2);
%     area18(:,:,i) = createMask(terr);
% end
% 
% save('../data/area18.mat', 'area18') % terrestrial area binary > 0mm rsl
% 
% treatment
% area19 = []; % initialize terrestrial area binary matrix
% for i= 1:nt_19 % loop through all timesteps
%     i
%     %What is sea level at time i
%     sl = (i*0.25*dt_19)+25; %first wet scan starts 48 minutes into the hour so i = i is a good approximation
%     elevationmask = ZW_19(:,:,i); 
%     elevationmask(elevationmask == 0) = NaN;
%     elevationmask_rslr = elevationmask - sl;
%     elevationmask_rslr(elevationmask_rslr < 0) = NaN;
% 
%     % Now make the shoreline boundary that removes any floating mats from
%     b = elevationmask_rslr; 
%     
%     % Use this for shoreline boundary
%     b(b > 0) = 1;
%     b(b <= 0) = 0;  
%     
%     b(isnan(b)) = 0;% address the nan values that become a problem later when creating a binary mask
%     B = bwboundaries(b,'noholes'); % Find the boundary of the matrix...creates n cell arrays of varying sizes
%     [msize, mindex] = max(cellfun('size',B,1)); % Find the cell array with the largest size...this is the cell array that contains the shoreline locations
%     C = B{mindex,1}; % Get the length of the largest cell array that contains the shoreline locations
% 
%     C2 = fliplr(C);
%     imagesc(b(:,:,1))
%     terr = drawpolygon('Position',C2);
%     area19(:,:,i) = createMask(terr);
% end
% 
% save('../data/area19.mat', 'area19', '-v7.3') % terrestrial area binary > 0mm rsl

%% Calculate the basin wide flow (total, channel, and overbank) properties 
% 1 is flow and 0 is no flow

% control 
flowscreen18 = flowscreen18.*area18; % multiply flow and terrestrial delta (1 is flow, 0 is no flow)

% initialize empty matricies
chanarea18 = [];
deltaarea18 = [];
flowarea18 = [];
obarea18 = [];
for i =1:nt_18
  chan = CM_18(:,:,i).*basin18;
  flow = flowscreen18(:,:,i).*basin18;
  flowtot = flow + chan; %deep parts of channel are sometimes not included in flow
  flowtot(flowtot>=1) = 1;
  ob = flowtot-chan; %0 channel or no flow, 1 ob flow, NaN outside basin or not terrestrial
  area = (area18(:,:,i)+chan).*basin18; % add channels since backwater length not captures in terrestrial area
  area(area>=1)=1; 
  chanarea18(i) = sum(chan(:), 'omitnan')*(2.5*10^-5); %m^2
  obarea18(i) = sum(ob(:), 'omitnan')*(2.5*10^-5); %m^2
  flowarea18(i) = sum(flowtot(:), 'omitnan')*(2.5*10^-5); %m^2
  deltaarea18(i) = sum(area(:), 'omitnan')*(2.5*10^-5); %m^2
end

% treatment 
% 1 is flow (and outside basin) and 0 is no flow
flowscreen19 = flowscreen19.*area19; % multiply flow, elevation relative to sea level, and terrestrial delta

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
  chan = CM_19(:,:,i).*basin19;
  area = (area19(:,:,i) + chan).*basin19; % add channels because backwater not included in terrestrial map
  area(area>=1)=1;
  % skip time steps with no channel map or cart impaired flow
  if sum(chan(:), 'omitnan') == 0 || ismember(i,idx) %if missing a channel map or cart impaired flow map, flow stats should be NaN
     chanarea19(:,i) = NaN;
     deltaarea19(:,i) = sum(area(:), 'omitnan')*(2.5*10^-5); %m2
     flowarea19(:,i) = NaN;
     obarea19(:,i) = NaN;
  else % calculate metrics for all other timesteps
      flow = flowscreen19(:,:,i).*basin19.*area;
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

%% Radial distance matrix
% Lets calculate the distance from apex to each pixel in our matrices
% control: create matrix of distances to apex
[X18, Y18] = meshgrid(1:ny_18, 1:nx_18); % x and y matrices
dd18 = sqrt((X18 - yentrance_18).^2 + (Y18 - xentrance_18).^2)*dx; % distance to pixel from apex; multiply by dx
% make everything outside of basin a NaN
tmp = zeros(nx_18,ny_18); % empty matrix
z = ZD_18(:,:,1); % intial z value to get basin boundaries
tmp2 = z.*tmp; % temporary matrix of basin; everything in basin = 0, everything outside = NaN
tmp2(tmp2 == 0.) = 1; % turn data inside basin into a 1
dd18 = dd18.*tmp2; % multiply data by distance to get a matrix with distance to eaach pixel (mm) for everything in the basin

% treatment: create matrix of distance to apex
[X19, Y19] = meshgrid(1:ny_19, 1:nx_19); % x and y matrices
dd19 = sqrt((X19 - yentrance_19).^2 + (Y19 - xentrance_19).^2)*dx; %distance to pixel from apex; multiply by dx (or dy)
%make everything outside of basin a NaN
tmp = zeros(nx_19,ny_19); % empty matrix to fill
z = ZD_19(:,:,1); % initial z value to get basin boundaries
z(z == 0.) = NaN; % remove well from matrix
tmp2 = z.*tmp; % temporary matrix of basin; everything in basin = 0; everything outside (and well) = NaN
tmp2(tmp2 == 0.) = 1; % turn data inside basin into a 1
dd19 = dd19.*tmp2; % multiply data by distance to get a matrix with distance to each pixel (mm) for everything in the basin (except the well)

%% Radial loop for channel and overbank area plots
rad_dist = 0:50:3100; % distance we will calculate for (from 0 to 3100 mm, by 50 mm or 5 cm).. we could change this if we want finer grained information

% initialize empty matrices 
chanfrac18_rad = NaN(length(rad_dist),nt_18); % radial channel fraction 
obfrac18_rad = NaN(length(rad_dist),nt_18); % number of channels

% loop through radial distances 
for k = 1:(length(rad_dist)-1) % loop to run through different radial distances from the apex.
    % print which radial transect the loop is on
    caption = sprintf('Radial segment %d', k);
    caption2 = sprintf(' of %d', length(rad_dist)-1);
    fprintf('%s\n', strcat(caption,caption2))
    
    % section to find x,y nodes for radial transect and generate matrix of
    % cross section topo and channel (yes/no) data
    radial_dd = dd18 >= rad_dist(k) & dd18 < rad_dist(k+1); % can look at this using imagesc(radial_dd) to visualize a 0.1 m radial transect

    % loop through time
    for i = 1:size(CM_18,3) % loop through all timesteps
        area = area18(:,:,i).*basin18.*radial_dd;
        area(area==0)=NaN;
        % channels
        chan = CM_18(:,:,i).*radial_dd;
        chan(chan == 0) = NaN;
        is_chan = (chan.*area); % in channel or no? on delta >0 mm or no?
        carea = sum(is_chan(:), 'omitnan')*0.25; % area of channel in cm^2
        
        % overbank
        is_chan(isnan(is_chan))=0;
        flow = flowscreen18(:,:,i).*radial_dd.*area; % in flow or no? on delta >0 mm or no?
        flowtot = flow + is_chan; %deep parts of channel are sometimes not included in flow
        flowtot(flowtot>=1) = 1;
        is_ob = flowtot-is_chan; %0 channel or no flow, 1 ob flow
        obarea = sum(is_ob(:), 'omitnan')*0.25; % area of overbank in cm^2

        % radial area
        rad = area.*radial_dd; % radial transect
        radial_area = sum(rad(:),'omitnan')*0.25; % area of radial transect in cm^2

        % save channel and overbank fraction
        chanfrac18_rad(k,i) = carea/radial_area;
        obfrac18_rad(k,i) = obarea/radial_area;

    end
end 

% treatment

% Visual spot checking shows we only need to remove the following 76 hours due to issues with cart messing up flow area:
idx = [79,85,87,89,91,93,95,97,99,...
    101,105,107,109,111,113,121,123,125,127,131,133,135,139,143,147,151,153,155,158,161,163,165,167,171,175,183,187,191,192,195,199,...
    203,207,215,219,227,230,231,239,243,247,248,251,255,259,267,271,291,293,294,297,...
    301,309,310,349,381,385,389,390,...
    401,405,417,443,449,450,...
    541];

% initialize empty matrices 
chanfrac19_rad = NaN(length(rad_dist),nt_19); % radial channel fraction 
obfrac19_rad = NaN(length(rad_dist),nt_19); % number of channels

% loop through radial distances 
for k = 1:(length(rad_dist)-1) % loop to run through different radial distances from the apex.
    % print which radial transect the loop is on
    caption = sprintf('Radial segment %d', k);
    caption2 = sprintf(' of %d', length(rad_dist)-1);
    fprintf('%s\n', strcat(caption,caption2))
    
    % section to find x,y nodes for radial transect and generate matrix of
    % cross section topo and channel (yes/no) data
    radial_dd = dd19 >= rad_dist(k) & dd19 < rad_dist(k+1); % can look at this using imagesc(radial_dd) to visualize a 0.1 m radial transect
    % skip time steps with no channel map or cart impaired flow

    % loop through time
    for i = 1:size(CM_19,3) % loop through all timesteps
        area = area19(:,:,i).*basin19.*radial_dd;
        area(area==0)=NaN;
        chan = CM_19(:,:,i);
        if sum(chan(:), 'omitnan') == 0 || ismember(i,idx) %if missing a channel map or cart impaired flow map, flow stats should be NaN
            chanfrac19_rad(:,i) = NaN;
            obfrac19_rad(:,i) = NaN;
        else 
            % channels
            chan = CM_19(:,:,i).*radial_dd;
            chan(chan == 0) = NaN;
            is_chan = (chan.*area); % in channel or no? on delta >0 mm or no?
            carea = sum(is_chan(:), 'omitnan')*0.25; % area of channel in cm^2
    
            % overbank
            is_chan(isnan(is_chan))=0;
            flow = flowscreen19(:,:,i).*radial_dd.*area; % in flow or no? on delta >0 mm or no?
            flowtot = flow + is_chan; %deep parts of channel are sometimes not included in flow
            flowtot(flowtot>=1) = 1;
            is_ob = flowtot-is_chan; %0 channel or no flow, 1 ob flow
            obarea = sum(is_ob(:), 'omitnan')*0.25; % area of overbank in cm^2
    
            % radial area
            rad = area.*radial_dd; % radial transect
            radial_area = sum(rad(:),'omitnan')*0.25; % area of radial transect in cm^2
    
            % save channel area and fraction
            chanfrac19_rad(k,i) = carea/radial_area;
            obfrac19_rad(k,i) = obarea/radial_area;   
        end
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
saveas(gcf,'../figures/esurf_Figure3a.pdf') % this was modified in illustrator for visual purposes

%% Plot the data (Figure 3b)
rad_dist = 0:0.05:3.1;
% create arrays
cf18 = mean(chanfrac18_rad,2, 'omitnan');
stdev18 = std(chanfrac18_rad,[],2, 'omitnan');
chan_array18 = [rad_dist; cf18'; stdev18'];
cols = any(isnan(chan_array18),1);
chan_array18(:,cols) = [];

cf19 = mean(chanfrac19_rad,2, 'omitnan');
stdev19 = std(chanfrac19_rad,[],2, 'omitnan');
chan_array19 = [rad_dist; cf19'; stdev19'];
cols = any(isnan(chan_array19),1);
chan_array19(:,cols) = [];

%fill standard deviation
y18 = chan_array18(2,:); % your mean vector;
x18 = chan_array18(1,:);
std18 = chan_array18(3,:);
curve1_18 = y18 + std18;
curve2_18 = y18 - std18;

y19 = chan_array19(2,:); % your mean vector;
x19 = chan_array19(1,:);
std19 = chan_array19(3,:);
curve1_19 = y19 + std19;
curve2_19 = y19 - std19;

of18 = mean(obfrac18_rad,2, 'omitnan');
stdev18 = std(obfrac18_rad,[],2, 'omitnan');
ob_array18 = [rad_dist; of18'; stdev18'];
cols = any(isnan(ob_array18),1);
ob_array18(:,cols) = [];

of19 = mean(obfrac19_rad,2, 'omitnan');
stdev19 = std(obfrac19_rad,[],2, 'omitnan');
ob_array19 = [rad_dist; of19'; stdev19'];
cols = any(isnan(ob_array19),1);
ob_array19(:,cols) = [];

%fill standard deviation
yob18 = ob_array18(2,:); % your mean vector;
xob18 = ob_array18(1,:);
std18 = ob_array18(3,:);
obcurve1_18 = yob18 + std18;
obcurve2_18 = yob18 - std18;

yob19 = ob_array19(2,:); % your mean vector;
xob19 = ob_array19(1,:);
std19 = ob_array19(3,:);
obcurve1_19 = yob19 + std19;
obcurve2_19 = yob19 - std19;

% Channel flow
fig = figure();
plot(x18, y18, 'b', 'LineWidth', 2)
hold on
plot(x19, y19, 'g', 'LineWidth', 2)
patch([x18 fliplr(x18)], [curve1_18 fliplr(curve2_18)], 'b')
patch([x19 fliplr(x19)], [curve1_19 fliplr(curve2_19)], 'g')
alpha(0.15)
plot(x18, y18, 'b', 'LineWidth', 2)
plot(x19, y19, 'g', 'LineWidth', 2)
ylim([0 0.8])
ylabel('channelized fraction (-)')
xlabel('distance from apex (m)') 
legend('control mean', 'treatment mean', 'control stdev', 'treatment stdev')
set(gca, 'XMinorTick', 'On', 'YMinorTick', 'On')
set(gcf, 'PaperUnits', 'inches');
y_width=7.25;x_width=9.125;
set(gcf, 'PaperPosition', [0 0 x_width y_width]);
saveas(fig, '../figures/esurf_Figure3b.pdf')

fig = figure;
plot(xob18, yob18, 'b', 'LineWidth', 2)
hold on
plot(xob19, yob19, 'g', 'LineWidth', 2)
patch([xob18 fliplr(xob18)], [obcurve1_18 fliplr(obcurve2_18)], 'b')
patch([xob19 fliplr(xob19)], [obcurve1_19 fliplr(obcurve2_19)], 'g')
alpha(0.15)
plot(xob18, yob18, 'b', 'LineWidth', 2)
plot(xob19, yob19, 'g', 'LineWidth', 2)
ylim([0 1])
ylabel('overbank fraction (-)')
xlabel('distance from apex (m)') 
legend('control mean', 'treatment mean', 'control stdev', 'treatment stdev')
set(gca, 'XMinorTick', 'On', 'YMinorTick', 'On')
set(gcf, 'PaperUnits', 'inches');
y_width=7.25;x_width=9.125;
set(gcf, 'PaperPosition', [0 0 x_width y_width]);
saveas(fig, '../figures/esurf_Figure3c.pdf')