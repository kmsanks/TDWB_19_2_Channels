clear all; close all
%% The purpose of this script is to calculate channel depth, width, and number of channels in the control and treatment experiments. 
% These data are plotted in Table 1, Figure 3a, 3c, and 5 in Sanks et al.
% (2023) submitted to Earth Surface Dynamics
% We will calculate the channel properties as a function of radial distance from the apex.

% Note that we will use every other hour of data for the control experiment
% in order to avoid any potential Saddler effects in our aggradation
% calculations, since the treatment experiment only contains LiDAR scans
% every other hour.

%% First, lets set our relative path
cd '../' %We should have opened this code from 'C:\Users\XXXXXXX\Documents\GitHub\TDWB_19_2_Channels\code'

%% Load data
cd './data'
load('ZD_18.mat') %control topography (mm)
load('ZD_19.mat') %treatment topography from Sam (mm)
load('ZD_19_2_dry.mat') %treatment topography from Jose (mm)
load('CM_18.mat') %control channel maps (binary, 1s = channels)
load('CM_19.mat') %treatment channel maps (binary, 1s = channels)
cd '../code'

%% Crop data to every two hours for both experiments to match ZD_19
% This step will allow us to calculate channel in-filling timescale
% since the control experiment has LiDAR every hour but the treatment is
% every other hour, we will change the control to match treat to avoid any
% potential saddler effects (specifically important for calcuating
% aggradation rates
ZD_18 = ZD_18(:,:,2:2:560);
CM_18 = CM_18(:,:,2:2:560);
%ZD_19 = ZD_19(:,:,2:end); % remove first time step since channel maps start at hour 2, and t = 1 is hour 0, t = 2 is hour 2
ZD_19 = ZD_19_2_dry(:,:,2:end); 
CM_19 = CM_19(:,:,2:2:end); % remove every other map, so we can start with hour 1
%% Define and set parameters
% control
% parameters
nx_18 = size(ZD_18, 1); %number of x locations on map
ny_18 = size(ZD_18,2); % number of y locations on map
nt_18 = size(ZD_18,3); % number of time steps in data set
dt_18 = 2; % delta t of time steps (hr)
xentrance_18 = 109; % x grid node location of the apex
yentrance_18 = 271; % y grid node location of the apex

% treatment
% parameters
nx_19 = size(ZD_19,1); % number of x locations on map
ny_19 = size(ZD_19,2); % number of y locations on map
nt_19 = size(ZD_19,3); % number of time steps in data set
dt_19 = 2; % delta t of time steps (hr)
xentrance_19 = 214; % x grid node location of the apex (x is down dip)
yentrance_19 = 397; % y grid node location of the apex (y is strike)

% both experiments
dx = 5; %5 mm grid cell in x
dy = 5; %5 mm grid cell in y
baselevel_rr = 0.25; % base level rise rate (mm/hr)
ocean_zero = 25; % ocean elevation at beginning of experiment (mm)

%% We only want to calculate for area above -9mm relative to sea level, so let's mask the basin
% control
terr_area18 = []; % initialize terrestrial area binary matrix
for i= 1:nt_18 % loop through all timesteps
    %What is sea level at time i
    sl = ((i+1)*0.25*dt_18)+25; % first scan starts at hour 2, so i+1
    elevationmask = ZD_18(:,:,i); 
    elevationmask(elevationmask == 0) = NaN;
    elevationmask_rslr = elevationmask - sl;
    elevationmask_rslr(elevationmask_rslr < -9) = NaN;

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

te

% treatment
terr_area19 = []; % initialize terrestrial area binary matrix
for i= 1:nt_19 % loop through all timesteps
    %What is sea level at time i
    sl = (i*0.25*dt_19)+25; % first scan starts at hour 1, so i
    elevationmask = ZD_19(:,:,i); 
    elevationmask(elevationmask == 0) = NaN;
    elevationmask_rslr = elevationmask - sl;
    elevationmask_rslr(elevationmask_rslr < -9) = NaN;

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

%% Extract channel topo
% Make elevation screen, so flow can be clipped to the areas above sea level
% control
tmp18 = [];
for i= 1:nt_18
    %What is sea level at time i
    sl = ((i+1)*baselevel_rr*dt_18)+25; %first wet scan starts at hour 2, so i+1
    elevationmask = ZD_18(:,:,i); 
    elevationmask(elevationmask == 0) = NaN;
    elevationmask_rslr = elevationmask - sl;
    elevationmask_rslr(elevationmask_rslr < -20) = NaN;
    tmp18(:,:,i) = elevationmask_rslr;
end

% treatment
tmp19 = [];
for i= 1:nt_19
    %What is sea level at time i
    sl = (i*baselevel_rr*dt_19)+25; %first scan is t = 2, which is hour 1, so i 
    elevationmask = ZD_19(:,:,i); 
    elevationmask(elevationmask == 0) = NaN;
    elevationmask_rslr = elevationmask - sl;
    elevationmask_rslr(elevationmask_rslr < -20) = NaN;
    tmp19(:,:,i) = elevationmask_rslr;
end

% Remove area outside channels 
% control
z18 = [];
for i =1:nt_18
  z18_RSL=(CM_18(:,:,i)).*(tmp18(:,:,i));
  z18_RSL(z18_RSL == 0) = NaN; %0 is outside channels 
  z18(:,:,i) = z18_RSL;
end

% treatment
z19 = [];
for i =1:nt_19
  z19_RSL=(CM_19(:,:,i)).*(tmp19(:,:,i));
  z19_RSL(z19_RSL == 0) = NaN; %0 is outside channels 
  z19(:,:,i) = z19_RSL;
end

%% Channel properties loop: here we will calculate radial channel properties (channel area, width, depth, number of channels)
% control first
dist = 0:50:3100; % distance we will calculate for (from 0 to 3100 mm, by 50 mm or 5 cm).. we could change this if we want finer grained information
% initialize empty matrices 
elev18 = NaN(length(dist),nt_18); % channel bed elevation relative to sea level
area18 = NaN(length(dist),nt_18); % total channel area
Hc_18 = NaN(length(dist),nt_18); % channel depth with distance and time to use in aggradation below
Hc_tot18 = NaN(length(dist),100); % to save percentile channel depth with distance
maxwidth18 = NaN(length(dist),nt_18); % trunk channel width
n_chan18 = NaN(length(dist),nt_18); % number of channels
% loop through radial distances 
for k = 1:(length(dist)-1) % loop to run through different radial distances from the apex.
    % section to find x,y nodes for radial transect and generate matrix of
    % cross section topo and channel (yes/no) data
    idx = dd18(dd18 >= dist(k) & dd18 < dist(k+1));
    radial_dd = dd18 >= dist(k) & dd18 < dist(k+1); % can look at this using imagesc(radial_dd) to visualize a 0.1 m radial transect
    % initialize empty matrices 
    depth = [];
    for i = 1:size(CM_18,3) % loop through all timesteps
        z = ZD_18(:,:,i); % elevation data for channel depth
        area = terr_area18(:,:,i) + CM_18(:,:,i);
        area(area>1)=1;
        is_shot = CM_18(:,:,i).*radial_dd.*area; % in channel or no? on delta >-9 mm or no?
        z = z.*is_shot; % elevation in the channels
        is_chan = is_shot;
        is_chan(is_chan == 0) = NaN;
        z_rsl = z18(:,:,i).*is_chan;
        bed = mean(z_rsl(:), 'omitnan');
        % save channel bed elevation for each distance through time
        if isempty(bed)
            elev18(k,i) = NaN;
        else 
            elev18(k,i) = bed;
        end
        [label,n] = bwlabel(is_shot); % label gives a unique number to each individual segement recognized; n is number of channels (or segments)
        ns_prop = regionprops(label, 'Area'); % this will calculate area for each of the segments
        list_area = (cell2mat({ns_prop.Area}'))*0.25; % area for each channel segment cm^2
        % save total channel area for each distance through time
        if isempty(list_area)
            area18(k,i) = NaN;
        else
            area18(k,i) = sum(list_area); 
        end
        % we will calculate depth for each channel here
        depthmax = [];
        for j = 1:n % loop through channel segments
            tmp = label;
            tmp(tmp ~= j) = NaN; % remove all data that is not in segement n
            tmp(tmp > 0) = 1; % turn segement data to 1
            tmp(tmp < 1) = NaN; % make everything else NaN
            z_tmp = z.*tmp; % segment elevations to get depth
            zs = z_tmp(~isnan(z_tmp)); % remove data not in the channel
            tmp_depth = (max(zs) - min(zs)); % channel depth in mm; we will take the max elevation as levee elevation and min elevation as thalweg                  
            depthmax = [depthmax, tmp_depth];
            depth = [depth,tmp_depth]; % to save all channel depths for each radial transect
        end
        depthmax = max(depthmax); % we don't want to average in the non-trunk channels, as they are often a lot shallower
        % save channel depth for each distance through time
        if isempty(depthmax)
            Hc_18(k,i) = NaN;
        else
            Hc_18(k,i) = depthmax;
        end
        list_width = (list_area)/5; % width in cm, assuming length is 5 mm
        list_width_max = max(list_width); % maximum width in radial segement (this is the "trunk" channel)
        % save trunk channel width for each distance through time
        if isempty(list_width_max)
            maxwidth18(k,i) = NaN;
        else
            maxwidth18(k,i) = list_width_max; 
        end
        % save number of channels for each distance through time
        if isempty(n)
            n_chan18(k,i) = NaN;
        else
            n_chan18(k,i) = n;
        end
    end
    for jj = 1:1:100
        Hc_tot18(k,jj) = prctile(depth,jj);
    end
end 

% treatment
% replace channel maps for treatment with NaN if time step does not have one
for i = 1:size(CM_19,3) 
    if sum(sum(CM_19(:,:,i), 'omitnan'), 'omitnan') == 0
        CM_19(:,:,i) = NaN;
    end
end

% initialize empty matrices
elev19 = NaN(length(dist),nt_19); % channel bed elevation relative to sea level
area19 = NaN(length(dist),nt_19); % total channel area
Hc_19 = NaN(length(dist),nt_19); % channel depth with distance and time to use in aggradation below
Hc_tot19 = NaN(length(dist),100); % to save 95th percentile channel depth with distance
maxwidth19 = NaN(length(dist),nt_19); % trunk channel width
n_chan19 = NaN(length(dist),nt_19); % number of channels
% run loop
for k = 1:(length(dist)-1)%loop to run through different radial distances from the end of the apex.
    %section to find x,y nodes for radial transect and generate matrix of
    %cross section topo and channel (yes/no) data
    idx = dd19(dd19 >= dist(k) & dd19 < dist(k+1));
    radial_dd = dd19 >= dist(k) & dd19 < dist(k+1);
    %how many channel pixels in the radius?
    depth = [];
    for i = 1:size(CM_19,3) % loop through each timestep
        z = ZD_19(:,:,i); % elevation data for channel depth
        area = terr_area19(:,:,i) + CM_19(:,:,i);
        area(area>1)=1;
        is_shot = CM_19(:,:,i).*radial_dd.*area; % in channel or no?
        z = z.*is_shot; % elevation in the channels
        if sum(sum(is_shot, 'omitnan')) == 0 % do not calculate channel properties for timesteps with no channel map
            elev19(k,i) = NaN;
            area19(k,i) = NaN;
            Hc_19(k,i) = NaN;
            maxwidth19(k,i) = NaN;
            n_chan19(k,i) = NaN;
        else %calculate channel properties for all other timesteps
            is_chan = is_shot;
            is_chan(is_chan == 0) = NaN;
            z_rsl = z19(:,:,i).*is_chan; 
            bed = mean(z_rsl(:), 'omitnan');
            % save channel bed elevation for each distance through time
            if isempty(bed)
                elev19(k,i) = NaN;
            else 
                elev19(k,i) = bed;
            end
            [label,n] = bwlabel(is_shot);
            ns_prop = regionprops(label, 'Area');
            list_area = (cell2mat({ns_prop.Area}'))*0.25; %cm^2
            % save total channel area for each distance through time
            if isempty(list_area)
                area19(k,i) = NaN;
            else
                area19(k,i) = sum(list_area); 
            end
            depthmax = [];
            for j = 1:n % loop through all channels
                tmp = label;
                tmp(tmp ~= j) = NaN; % remove all data that is not in segement n
                tmp(tmp > 0) = 1; % turn segement data to 1
                tmp(tmp < 1) = NaN; % make everything else NaN
                z_tmp = z.*tmp; % segment elevations to get depth
                zs = z_tmp(~isnan(z_tmp)); % remove data not in the channel
                tmp_depth = (max(zs) - min(zs)); % channel depth in mm; we will take the max elevation as levee elevation and min elevation as thalweg         
                depthmax = [depthmax, tmp_depth];
                depth = [depth,tmp_depth]; % to save all channel depths for each radial transect
            end
            depthmax = max(depthmax); % we don't want to average in the non-trunk channels, as they are often a lot shallower
            if isempty(depthmax)
                Hc_19(k,i) = NaN;
            else
                Hc_19(k,i) = depthmax;
            end
            list_width = (list_area)/5; % divide by channel length, which is 5 cm
            list_width_max = max(list_width);
            % save trunk channel width for each distance through time
            if isempty(list_width_max)
                maxwidth19(k,i) = NaN;
            else
                maxwidth19(k,i) = list_width_max; 
            end
            % save number of channels for each distance through time
            if isempty(n)
                n_chan19(k,i) = NaN;
            else
                n_chan19(k,i) = n;
            end
        end
    end
    for jj = 1:1:100
        Hc_tot19(k,jj) = prctile(depth,jj);
    end
end 

%% Calculate statistics for Table 1:
% trunk channel width
width_mean18 = mean(maxwidth18(:), 'omitnan');
width_std18 = std(maxwidth18(:), 'omitnan');
width_mean19 = mean(maxwidth19(:), 'omitnan');
width_std19 = std(maxwidth19(:), 'omitnan');

% trunk channel depth
depth_mean18 = mean(Hc_18(:), 'omitnan');
depth_std18 = std(Hc_18(:), 'omitnan');
depth_mean19 = mean(Hc_19(:), 'omitnan');
depth_std19 = std(Hc_19(:), 'omitnan');

% trunk channel depth without entrance
tmpdepth = Hc_18(6:63,:);
depth_mean18_noent = mean(tmpdepth(:), 'omitnan');
depth_std18_noent = std(tmpdepth(:), 'omitnan');
tmpdepth = Hc_19(6:63,:);
depth_mean19_noent = mean(tmpdepth(:), 'omitnan');
depth_std19_noent = std(tmpdepth(:), 'omitnan');

% mean channel bed elevation relative to sea level
celev_mean18 = mean(elev18(:), 'omitnan');
celev_std18 = std(elev18(:), 'omitnan');
celev_mean19 = mean(elev19(:), 'omitnan');
celev_std19 = std(elev19(:), 'omitnan');

%% Plot the data: Figure 3a
ybars = [-9 5]; % marsh window
dist = 0:0.05:3.1;
el18 = mean(elev18,2, 'omitnan');
stdev18 = std(elev18,[],2, 'omitnan');
array18 = [dist; el18'; stdev18'];
cols = any(isnan(array18),1);
array18(:,cols) = [];

el19 = mean(elev19,2, 'omitnan');
stdev19 = std(elev19,[],2, 'omitnan');
array19 = [dist; el19'; stdev19'];
cols = any(isnan(array19),1);
array19(:,cols) = [];

%fill standard deviation
y18 = array18(2,:); % your mean vector;
x18 = array18(1,:);
std18 = array18(3,:);
curve1_18 = y18 + std18;
curve2_18 = y18 - std18;

y19 = array19(2,:); % your mean vector;
x19 = array19(1,:);
std19 = array19(3,:);
curve1_19 = y19 + std19;
curve2_19 = y19 - std19;

fig = figure();
plot(x18, y18, 'b', 'LineWidth', 2)
hold on
plot(x19, y19, 'g', 'LineWidth', 2)
patch([x18 fliplr(x18)], [curve1_18 fliplr(curve2_18)], 'b')
patch([x19 fliplr(x19)], [curve1_19 fliplr(curve2_19)], 'g')
patch([min(xlim) max(xlim) max(xlim) min(xlim)], [ybars(1) ybars(1), ybars(2) ybars(2)], 'k')
alpha(0.15)
yline(0, 'k-', 'linewidth', 2)
plot(x18, y18, 'b', 'LineWidth', 2)
plot(x19, y19, 'g', 'LineWidth', 2)
ylim([-10 40])
xlim([0 3])
ylabel('mean channel bed elevation relative to sea level (mm)')
xlabel('distance from apex (m)') 
legend('control mean', 'treatment mean', 'control stdev', 'treatment stdev', 'marsh window', 'sea level')
set(gca, 'XMinorTick', 'On', 'YMinorTick', 'On')
set(gcf, 'PaperUnits', 'inches');
y_width=7.25;x_width=9.125;
set(gcf, 'PaperPosition', [0 0 x_width y_width]);
saveas(fig, '../figures/esurf_Figure3a.pdf')

%% Channel depth: visualization
d18 = mean(Hc_18,2, 'omitnan');
stdev18 = std(Hc_18,[],2, 'omitnan');
array18 = [dist; d18'; stdev18'];
cols = any(isnan(array18),1);
array18(:,cols) = [];

d19 = mean(Hc_19,2, 'omitnan');
stdev19 = std(Hc_19,[],2, 'omitnan');
array19 = [dist; d19'; stdev19'];
cols = any(isnan(array19),1);
array19(:,cols) = [];

%fill standard deviation
y18 = array18(2,:); % your mean vector;
x18 = array18(1,:);
std18 = array18(3,:);
curve1_18 = y18 + std18;
curve2_18 = y18 - std18;

y19 = array19(2,:); % your mean vector;
x19 = array19(1,:);
std19 = array19(3,:);
curve1_19 = y19 + std19;
curve2_19 = y19 - std19;

figure()
plot(x18, y18, 'b', 'LineWidth', 2)
hold on
plot(x19, y19, 'g', 'LineWidth', 2)
patch([x18 fliplr(x18)], [curve1_18 fliplr(curve2_18)], 'b')
patch([x19 fliplr(x19)], [curve1_19 fliplr(curve2_19)], 'g')
alpha(0.15)
plot(x18, y18, 'b', 'LineWidth', 2)
plot(x19, y19, 'g', 'LineWidth', 2)
xlim([0 3])
ylabel('channel depth (mm)')
xlabel('radial distance from channel entrance (m)') 
legend('control mean', 'treatment mean', 'control stdev', 'treatment stdev')
%set(gcf, 'PaperUnits', 'inches');
%y_width=7.25;x_width=9.125;
%set(gcf, 'PaperPosition', [0 0 x_width y_width]);
%saveas(fig, '../figures/esurf_Figure3a.pdf')


%% Plot the data: Figure 3c
% create arrays
width18 = mean(maxwidth18,2, 'omitnan');
stdev18 = std(maxwidth18,[],2, 'omitnan');
chan_array18 = [dist; width18'; stdev18'];
cols = any(isnan(chan_array18),1);
chan_array18(:,cols) = [];

width19 = mean(maxwidth19,2, 'omitnan');
stdev19 = std(maxwidth19,[],2, 'omitnan');
chan_array19 = [dist; width19'; stdev19'];
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

% make figure
fig = figure();
yyaxis left
patch([x18 fliplr(x18)], [curve1_18 fliplr(curve2_18)], 'b')
hold on
patch([x19 fliplr(x19)], [curve1_19 fliplr(curve2_19)], 'g')
alpha(0.15)
plot(x18, y18, 'b-', 'LineWidth', 2)
plot(x19, y19, 'g-', 'LineWidth', 2)
ylabel('max channel width (cm)')
ylim([0 18])
xlim([0 3])
set(gca, 'YMinorTick', 'On')
yyaxis right
plot(0:0.05:3.1, mean(n_chan18,2,'omitnan'), 'bx')
plot(0:0.05:3.1, mean(n_chan19,2,'omitnan'), 'gx')
ylim([0 4])
xlim([0 3])
ylabel('number of channels')
xlabel('distance from apex (m)') 
legend('control mean', 'treatment mean', 'control stdev', 'treatment stdev', 'number of channels control', 'number of channels treatment')
ax = gca;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'k';
set(gca, 'YMinorTick', 'On', 'XMinorTick', 'On', 'XAxisLocation', 'bottom', 'XAxisLocation', 'bottom')
set(gcf, 'PaperUnits', 'inches');
y_width=7.25;x_width=9.125;
set(gcf, 'PaperPosition', [0 0 x_width y_width]);
saveas(fig, '../figures/esurf_Figure3c.pdf')

%% Now we will calculate channel in-filling and compensation
% reload treatment channel maps?

% Process the elevation data to remove extranneous values
% control
nt_18 = size(CM_18, 3);
z18 = []; %go by 2 here so no saddler effects when comparing 18 and 19
for i = 1:nt_18
    z = ZD_18(:,:,i);
    z_rsl = z - (baselevel_rr*i*dt_18+ocean_zero);
    z_rsl(z_rsl<0) = NaN;
    z18(:,:,i) = z_rsl + (baselevel_rr*i*dt_18+ocean_zero);
end 

% treatment
z19 = [];
for i = 1:nt_19
    z = ZD_19(:,:,i);
    z_rsl = z - (baselevel_rr*(i-1)*dt_19+ocean_zero); %subtract 1 from i since hour 1 here is technically hour 2
    z_rsl(z_rsl<0) = NaN;
    z19(:,:,i) = z_rsl + (baselevel_rr*(i-1)*dt_19+ocean_zero);
end 

%% Lets buffer the channels so levee sedimentation is included in channel aggradation 
se = strel('square',10); %create a square of 10 pixels around each 1 (channel) pixel, this will buffer our channel maps
% control
CM_buffer_18 = [];
for j = 1:nt_18 %go by 2 so no saddler effects here
    buffer = imdilate(CM_18(:,:,j), se); %buffer channel map by 10 pixels
    CM_buffer_18(:,:,j) = buffer;
end 
% treatment
CM_buffer_19 = [];
for j = 1:nt_19
    buffer = imdilate(CM_19(:,:,j), se);
    CM_buffer_19(:,:,j) = buffer;
end 

%% Clean treatment data
% Replace timesteps with no channel maps with the channel map from the
% next time step for the treatment experiment
for i = (nt_19-1):-1:1 %I know 280 has a channel map, so I can start at 279 and replace with channel map that comes after, this will work for ones that have multiple no maps in a row
    if sum(sum(CM_buffer_19(:,:,i), 'omitnan'), 'omitnan') == -Inf
        CM_buffer_19(:,:,i) = CM_buffer_19(:,:,i+1);
    end
end

%% Calculate channel and far-field aggradation
% Total dz: control
dZD_18 = diff(ZD_18,1,3); % difference the maps along time dimension (3)
% Total dz: treatment
dZD_19 = diff(ZD_19,1,3); % difference the maps along time dimension (3)

% Far field maps (remove channel from topo): control
FF_18_buffer = ~CM_buffer_18; % inverse of channel buffer
ZD_18_FF_buffer = FF_18_buffer.*z18; % elevation of far field
ZD_18_FF_buffer(ZD_18_FF_buffer == 0) = NaN; % turn channel area into NaN
dZD_18_FF_buffer = diff(ZD_18_FF_buffer,1,3); % difference the maps along time dimension (3)

% Far field maps (remove channel from topo): treatment
FF_19_buffer = ~CM_buffer_19; % inverse of channel buffer
ZD_19_FF_buffer = FF_19_buffer.*z19; % elevation of far field
ZD_19_FF_buffer(ZD_19_FF_buffer == 0) = NaN; % turn channel area into NaN
dZD_19_FF_buffer = diff(ZD_19_FF_buffer,1,3); % difference the maps along time dimension (3)

% Channel maps (remove far field from topo): control
ZD_18_chan_buffer = CM_buffer_18.*z18; % elevation of channels
ZD_18_chan_buffer(ZD_18_chan_buffer == 0) = NaN; % remove far field area 
dZD_18_chan_buffer = diff(ZD_18_chan_buffer,1,3); % difference the maps along time dimension (3)

% Channel maps (remove far field from topo): treatment
ZD_19_chan_buffer = CM_buffer_19.*z19; % elevation of channels
ZD_19_chan_buffer(ZD_19_chan_buffer == 0) = NaN; % remove far field area
dZD_19_chan_buffer = diff(ZD_19_chan_buffer,1,3); % difference the maps along time dimension (3)

%% Calculate aggradation statistics for Table 1
mean_chan_agg18 = mean(dZD_18_chan_buffer(:), 'omitnan')/2; %mm/hr
stdev_chan_agg18 = std(dZD_18_chan_buffer(:), 'omitnan')/2; %mm/hr

mean_chan_agg19 = mean(dZD_19_chan_buffer(:), 'omitnan')/2; %mm/hr
stdev_chan_agg19 = std(dZD_19_chan_buffer(:), 'omitnan')/2; %mm/hr

mean_ff_agg18 = mean(dZD_18_FF_buffer(:), 'omitnan')/2; %mm/hr
stdev_ff_agg18 = std(dZD_18_FF_buffer(:), 'omitnan')/2; %mm/hr

mean_ff_agg19 = mean(dZD_19_FF_buffer(:), 'omitnan')/2; %mm/hr
stdev_ff_agg19 = std(dZD_19_FF_buffer(:), 'omitnan')/2; %mm/hr

%% Calculate aggradation, compensation timescale, and channel in-filling rate
% control
CM_buffer_18 = CM_buffer_18(:,:,1:279); % when you difference the maps, you will have one less time step 
FF_18_buffer = FF_18_buffer(:,:,1:279); % far field maps

% initialize empty matrices to fill 
dist = 0:50:3100; % each radial transect is 50 mm or 5 cm
agg18 = NaN(length(dist),nt_18-1); % total aggradation rate for distance through time
c_agg18 = NaN(length(dist),nt_18-1); % channel aggradation rate for distance through time
ff_agg18 = NaN(length(dist),nt_18-1); % far-field aggradation rate for distance through time
d_agg18 = NaN(length(dist),nt_18-1); % difference in mean aggradation between channel and ff for distance through time
% loop through radial distances 
for k = 1:(length(dist)-1) % loop to run through different radial distances from the end of the entrance channel.
    % section to find x,y nodes for radial transect and generate matrix of
    % cross section topo and channel (yes/no) data
    idx = dd18(dd18 >= dist(k) & dd18 < dist(k+1));
    radial_dd = dd18 >= dist(k) & dd18 < dist(k+1); % can look at this using imagesc(radial_dd) to visualize a 0.1 m radial transect
    % initialize empty matrices 
    for i = 1:(size(CM_18,3)-1)
        % only want area on the delta top >-9 mm rsl or in the channels 
        area = terr_area18(:,:,i) + CM_18(:,:,i);
        area(area>1)=1;
        % dz
        dz = dZD_18(:,:,i).*radial_dd.*area;
        dz(dz == 0.) = NaN; % remove data outside radial transect
        % channel and ff data in radial transect
        dz_chan = dZD_18_chan_buffer(:,:,i).*radial_dd.*area; % multiply by radial transect
        dz_chan(dz_chan == 0.) = NaN; %remove channel pixels not in radial_dd
        dz_ff = dZD_18_FF_buffer(:,:,i).*radial_dd.*area; % multiply by radial transect
        dz_ff(dz_ff == 0.) = NaN; % remove ff pixels not in radial_dd
        % now we will calculate aggradation rates
        % total aggradation
        dz_t = dz(~isnan(dz(:)));
        if isempty(dz_t)
            agg18(k,i) = NaN;
        else
            agg18(k,i) = mean(dz_t/2, 'omitnan'); % mm/hr
        end
        % channel aggradation
        dz_c = dz_chan(~isnan(dz_chan(:)));
        if isempty(dz_c)
            c_agg18(k,i) = NaN;
        else 
            c_agg18(k,i) = mean(dz_c/2, 'omitnan'); % mm/hr
        end
        % far-field aggradation
        dz_ff = dz_ff(~isnan(dz_ff(:)));
        if isempty(dz_ff)
            ff_agg18(k,i) = NaN;
        else
            ff_agg18(k,i) = mean(dz_ff/2, 'omitnan'); % mm/hr
        end
        % difference in agg between channel and ff
        if isempty(dz_c) || isempty(dz_ff)
            d_agg18(k,i) = NaN;
        else
            d_agg18(k,i) = mean(dz_c/2, 'omitnan') - mean(dz_ff/2, 'omitnan'); % mm/hr
        end
    end
end 

% treatment
CM_buffer_19 = CM_buffer_19(:,:,1:279); % when you difference the maps, you will have one less time step 
FF_19_buffer = FF_19_buffer(:,:,1:279); % far field maps

% initialize empty matrices to fill  
agg19 = NaN(length(dist),nt_19-1); % total aggradation rate for distance through time
c_agg19 = NaN(length(dist),nt_19-1); % channel aggradation rate for distance through time
ff_agg19 = NaN(length(dist),nt_19-1); % far-field aggradation rate for distance through time
d_agg19 = NaN(length(dist),nt_19-1); % difference in mean aggradation between channel and ff for distance through time
% loop through radial distances 
for k = 1:(length(dist)-1) % loop to run through different radial distances from the end of the entrance channel.
    % section to find x,y nodes for radial transect and generate matrix of
    % cross section topo and channel (yes/no) data
    idx = dd19(dd19 >= dist(k) & dd19 < dist(k+1));
    radial_dd = dd19 >= dist(k) & dd19 < dist(k+1); % can look at this using imagesc(radial_dd) to visualize a 0.1 m radial transect
    for i = 1:(size(CM_19,3)-1)
        % only want area on the delta top >-9 mm rsl or in the channels 
        area = terr_area19(:,:,i) + CM_19(:,:,i);
        area(area>1)=1;
        % dz
        dz = dZD_19(:,:,i).*radial_dd.*area;
        dz(dz == 0.) = NaN; % remove data outside radial transect
        % channel and ff data in radial transect
        dz_chan = dZD_19_chan_buffer(:,:,i).*radial_dd.*area; % multiply by radial transect
        dz_chan(dz_chan == 0.) = NaN; %remove channel pixels not in radial_dd
        dz_ff = dZD_19_FF_buffer(:,:,i).*radial_dd.*area; % multiply by radial transect
        dz_ff(dz_ff == 0.) = NaN; % remove ff pixels not in radial_dd
        % now we will calculate aggradation rates
        % total aggradation
        dz_t = dz(~isnan(dz(:)));
        if isempty(dz_t)
            agg19(k,i) = NaN;
        else
            agg19(k,i) = mean(dz_t/2, 'omitnan'); % mm/hr
        end
        % channel aggradation
        dz_c = dz_chan(~isnan(dz_chan(:)));
        if isempty(dz_c)
            c_agg19(k,i) = NaN;
        else 
            c_agg19(k,i) = mean(dz_c/2, 'omitnan'); % mm/hr
        end
        % far-field aggradation
        dz_ff = dz_ff(~isnan(dz_ff(:)));
        if isempty(dz_ff)
            ff_agg19(k,i) = NaN;
        else
            ff_agg19(k,i) = mean(dz_ff/2, 'omitnan'); % mm/hr
        end
        % difference in agg between channel and ff
        if isempty(dz_c) || isempty(dz_ff)
            d_agg19(k,i) = NaN;
        else
            d_agg19(k,i) = mean(dz_c/2, 'omitnan') - mean(dz_ff/2, 'omitnan'); % mm/hr
        end
    end 
end 

%% Calculate compensation timescale and channel in-filling rate
hc = Hc_18(:,1:279);
hc18 = mean(hc,2, 'omitnan');
agg18_tmp = mean(agg18,2, 'omitnan');
dagg18 = mean(d_agg18,2, 'omitnan');
Tc_18 = hc18./agg18_tmp;
Ta_18 = hc18./dagg18;

hc = Hc_19(:,1:279);
hc19 = mean(hc,2, 'omitnan');
agg19_tmp = mean(agg19,2, 'omitnan');
dagg19 = mean(d_agg19,2, 'omitnan');
Tc_19 = hc19./agg19_tmp;
Ta_19 = hc19./dagg19;
%infill_std = infill*(sqrt((dagg_std/dagg)^2+(depth_std18(k)/depth18(k))^2)); % error propagation
% control
% dagg_std = sqrt((std(c_agg18)^2)+(std(f_agg18)^2)/2); % error propagation


%% Plot the data: compensation timescale
% create arrays
dist = 0:0.05:3.1;
tc18 = mean(Tc_18,2, 'omitnan');
stdev18 = std(Tc_18,[],2, 'omitnan');
chan_array18 = [dist; tc18'; stdev18'];
cols = any(isnan(chan_array18),1);
chan_array18(:,cols) = [];

tc19 = mean(Tc_19,2, 'omitnan');
stdev19 = std(Tc_19,[],2, 'omitnan');
chan_array19 = [dist; tc19'; stdev19'];
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

% make figure
fig = figure();
patch([x18 fliplr(x18)], [curve1_18 fliplr(curve2_18)], 'b')
hold on
patch([x19 fliplr(x19)], [curve1_19 fliplr(curve2_19)], 'g')
alpha(0.15)
plot(x18, y18, 'b-', 'LineWidth', 2)
plot(x19, y19, 'g-', 'LineWidth', 2)
ylabel('compensation timescale (hrs)')
ylim([0 300])
legend('control mean', 'treatment mean', 'control stdev', 'treatment stdev')
set(gca, 'XMinorTick', 'On', 'YMinorTick', 'On')
set(gcf, 'PaperUnits', 'inches');
y_width=7.25;x_width=9.125;
set(gcf, 'PaperPosition', [0 0 x_width y_width]);
saveas(fig, '../figures/esurf_Figure5X.pdf')
%% Plot the data: avulsion timescale
% create arrays
dist = 0:0.05:3.1;
ta18 = mean(Ta_18,2, 'omitnan');
stdev18 = std(Ta_18,[],2, 'omitnan');
chan_array18 = [dist; ta18'; stdev18'];
cols = any(isnan(chan_array18),1);
chan_array18(:,cols) = [];

ta19 = mean(Ta_19,2, 'omitnan');
stdev19 = std(Ta_19,[],2, 'omitnan');
chan_array19 = [dist; ta19'; stdev19'];
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

fig = figure();
plot(x18, y18, 'b-', 'LineWidth', 2)
hold on
plot(x19, y19, 'g-', 'LineWidth', 2)
ylabel('avulsion timescale (hrs)')
ylim([0 300])
legend('control mean', 'treatment mean')
set(gca, 'XMinorTick', 'On', 'YMinorTick', 'On')
set(gcf, 'PaperUnits', 'inches');
y_width=7.25;x_width=9.125;
set(gcf, 'PaperPosition', [0 0 x_width y_width]);
saveas(fig, '../figures/esurf_Figure5c.pdf')
%% Plot the data: Figure 5
% error bars for plotting
% channel aggradation
chan_agg_mean_18 = mean(c_agg18,2, 'omitnan');
chan_agg_std_18 = mean(c_agg18,2, 'omitnan');
chan_array18 = [0:0.05:3.1; chan_agg_mean_18'; chan_agg_std_18'];
cols = any(isnan(chan_array18),1);
chan_array18(:,cols) = [];

chan_agg_mean_19 = mean(c_agg19,2, 'omitnan');
chan_agg_std_19 = mean(c_agg19,2, 'omitnan');
chan_array19 = [0:0.05:3.1; chan_agg_mean_19'; chan_agg_std_19'];
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

% ff aggradation
ff_agg_mean_18 = mean(ff_agg18,2, 'omitnan');
ff_agg_std_18 = mean(ff_agg18,2, 'omitnan');
ff_array18 = [0:0.05:3.1; ff_agg_mean_18'; ff_agg_std_18'];
cols = any(isnan(ff_array18),1);
ff_array18(:,cols) = [];

ff_agg_mean_19 = mean(ff_agg19,2, 'omitnan');
ff_agg_std_19 = mean(ff_agg19,2, 'omitnan');
ff_array19 = [0:0.05:3.1; ff_agg_mean_19'; ff_agg_std_19'];
cols = any(isnan(ff_array19),1);
ff_array19(:,cols) = [];

%fill standard deviation
ff_y18 = ff_array18(2,:); % your mean vector;
ff_x18 = ff_array18(1,:);
ff_std18 = ff_array18(3,:);
ff_curve1_18 = ff_y18 + ff_std18;
ff_curve2_18 = ff_y18 - ff_std18;

ff_y19 = ff_array19(2,:); % your mean vector;
ff_x19 = ff_array19(1,:);
ff_std19 = ff_array19(3,:);
ff_curve1_19 = ff_y19 + ff_std19;
ff_curve2_19 = ff_y19 - ff_std19;

% difference in agg
% dd aggradation
d_agg_mean_18 = mean(dagg18,2, 'omitnan');
d_agg_std_18 = mean(dagg18,2, 'omitnan');
dd_array18 = [0:0.05:3.1; d_agg_mean_18'; d_agg_std_18'];
cols = any(isnan(dd_array18),1);
dd_array18(:,cols) = [];

d_agg_mean_19 = mean(dagg19,2, 'omitnan');
d_agg_std_19 = mean(dagg19,2, 'omitnan');
dd_array19 = [0:0.05:3.1; d_agg_mean_19'; d_agg_std_19'];
cols = any(isnan(dd_array19),1);
dd_array19(:,cols) = [];

%fill standard deviation
dd_y18 = dd_array18(2,:); % your mean vector;
dd_x18 = dd_array18(1,:);
dd_std18 = dd_array18(3,:);
dd_curve1_18 = dd_y18 + dd_std18;
dd_curve2_18 = dd_y18 - dd_std18;

dd_y19 = dd_array19(2,:); % your mean vector;
dd_x19 = dd_array19(1,:);
dd_std19 = dd_array19(3,:);
dd_curve1_19 = dd_y19 + dd_std19;
dd_curve2_19 = dd_y19 - dd_std19;

% Figure 5a
% channel agg
fig = figure();
plot(x18, y18, 'b-', 'LineWidth', 2)
hold on
plot(x19, y19, 'g-', 'LineWidth', 2)
plot(ff_x18, ff_y18, 'b--', 'LineWidth', 2)
plot(ff_x19, ff_y19, 'g--', 'LineWidth', 2)
yline(0.25,'k:', 'LineWidth', 2)
ylabel('aggradation rate (mm/hr)')
xlabel('distance from apex (m)')
legend('control channel', 'treatment channel', 'control far-field', 'treatment far-field', 'RSLR{_b}')
set(gca,'XMinorTick','on','YMinorTick','on')
y_width=7.25;x_width=9.125;
set(gcf, 'PaperPosition', [0 0 x_width y_width]);
saveas(fig, '../figures/esurf_Figure5a.pdf')

% Figure 5b
fig = figure();
patch([dd_x18 fliplr(dd_x18)], [dd_curve1_18 fliplr(dd_curve2_18)], 'b--')
hold on
patch([dd_x19 fliplr(dd_x19)], [dd_curve1_19 fliplr(dd_curve2_19)], 'g--')
alpha(0.15)
plot(dd_x18, dd_y18, 'b-', 'LineWidth', 2)
plot(dd_x19, dd_y19, 'g-', 'LineWidth', 2)
ylabel('channel in-filling rate (mm/hr)')
ylim([0 3])
xlabel('distance from apex (m)') 
legend('control mean', 'treatment mean', 'control stdev', 'treatment stdev') 
set(gca,'XMinorTick','on','YMinorTick','on')
y_width=7.25;x_width=9.125;
set(gcf, 'PaperPosition', [0 0 x_width y_width]);
saveas(fig, '../figures/esurf_Figure5b.pdf')