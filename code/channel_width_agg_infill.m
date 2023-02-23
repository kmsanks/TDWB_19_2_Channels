clear all; close all
%% The purpose of this script is to calculate channel depth, width, and number of channels in the control and treatment experiments. 
% These data are plotted in Table 1 and Figure 3c in Sanks et al.
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
load('ZD_19.mat') %treatment topography (mm)
load('CM_18.mat') %control channel maps (binary, 1s = channels)
load('CM_19.mat') %treatment channel maps (binary, 1s = channels)
cd '../code'

%% Define and set parameters
% control
% parameters
nx_18 = size(ZD_18, 1); %number of x locations on map
ny_18 = size(ZD_18,2); % number of y locations on map
nt_18 = size(ZD_18,3); % number of time steps in data set
dt_18 = 1; % delta t of time steps (hr)
xentrance_18 = 109; % x grid node location of the apex
yentrance_18 = 271; % y grid node location of the apex

% treatment
CM_19 = CM_19(:,:,2:2:end); % change channel maps to match LiDAR
ZD_19 = ZD_19(:,:,2:end); % remove first time step since channel maps start at hour 2, and t = 1 is hour 0, t = 2 is hour 2
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

%% Calculate channel width
% Lets calculate the distance from apex to each pixel in our matrices
% control: create matrix of distances to apex
[X18 Y18] = meshgrid(1:ny_18, 1:nx_18); % x and y matrices
dd18 = sqrt((X18 - yentrance_18).^2 + (Y18 - xentrance_18).^2)*dx; % distance to pixel from apex; multiply by dx
% make everything outside of basin a NaN
tmp = zeros(nx_18,ny_18); % empty matrix
z = ZD_18(:,:,1); % intial z value to get basin boundaries
tmp2 = z.*tmp; % temporary matrix of basin; everything in basin = 0, everything outside = NaN
tmp2(tmp2 == 0.) = 1; % turn data inside basin into a 1
dd18 = dd18.*tmp2; % multiply data by distance to get a matrix with distance to eaach pixel (mm) for everything in the basin

% treatment: create matrix of distance to apex
[X19 Y19] = meshgrid(1:ny_19, 1:nx_19); % x and y matrices
dd19 = sqrt((X19 - yentrance_19).^2 + (Y19 - xentrance_19).^2)*dx; %distance to pixel from apex; multiply by dx (or dy)
%make everything outside of basin a NaN
tmp = zeros(nx_19,ny_19); % empty matrix to fill
z = ZD_19(:,:,1); % initial z value to get basin boundaries
z(z == 0.) = NaN; % remove well from matrix
tmp2 = z.*tmp; % temporary matrix of basin; everything in basin = 0; everything outside (and well) = NaN
tmp2(tmp2 == 0.) = 1; % turn data inside basin into a 1
dd19 = dd19.*tmp2; % multiply data by distance to get a matrix with distance to each pixel (mm) for everything in the basin (except the well)

%% Channel properties loop: here we will calculate radial channel properties (channel area, width, depth, number of channels)
% control first
% initialize empty matrices 
area18 = []; % channel area
area_std18 =[]; % standard deviation of channel area
depth18 = []; % channel depth
depth_std18 = []; % standard deviation of channel depth 
width18 = []; % channel width
width_std18 = []; % standard deviation of channel width
maxwidth18 = []; % maximum channel width
maxwidth_std18 = []; % standard deviation of maximum channel width
n_chan18 = []; % number of channels
n_chan_std18 =[]; % standard deviation of the number of channels
dist = 0:100:3100; % distance we will calculate for (from 0 to 3100 mm, by 100 mm or 0.1 m).. we could change this if we want finer grained information

% loop through radial distances 
for k = 1:(length(dist)-1) % loop to run through different radial distances from the apex.
    % section to find x,y nodes for radial transect and generate matrix of
    % cross section topo and channel (yes/no) data
    idx = dd18(dd18 >= dist(k) & dd18 < dist(k+1));
    radial_dd = dd18 >= dist(k) & dd18 < dist(k+1); % can look at this using imagesc(radial_dd) to visualize a 0.1 m radial transect
    % initialize empty matrices 
    chan_area18 = []; 
    chan_depth18 = [];
    chan_width18 = [];
    chan_width_max18 = [];
    n_chan_i = [];
    for i = 1:size(CM_18,3) % loop through all timesteps
        z = ZD_18(:,:,i); % elevation data for channel depth
        is_shot = CM_18(:,:,i).*radial_dd; % in channel or no?
        z = z.*is_shot; % elevation in the channels
        [label,n] = bwlabel(is_shot); % label gives a unique number to each individual segement recognized; n is number of channels (or segments)
        ns_prop = regionprops(label, 'Area'); % this will calculate area for each of the segments
        list_area = (cell2mat({ns_prop.Area}'))*0.25; % area for each channel segment cm^2
        % we will calculate depth for each channel here
        depth = [];
        for j = 1:n % loop through channel segments
            tmp = label;
            tmp(tmp ~= j) = NaN; % remove all data that is not in segement n
            tmp(tmp > 0) = 1; % turn segement data to 1
            tmp(tmp < 1) = NaN; % make everything else NaN
            z_tmp = z.*tmp; % segment elevations to get depth
            zs = z_tmp(~isnan(z_tmp)); % remove data not in the channel
            tmp_depth = (prctile(zs, 95) - prctile(zs, 5)); % channel depth in mm; we will take the 5% value as our channel bottom and the 95% value as our levee elevation, subtract the two.         
            depth = [depth, tmp_depth];
        end
        depth = max(depth); % we don't want to average in the non-trunk channels, as they are often a lot shallower
        chan_depth18 = [chan_depth18;depth];
        %list_diam = cell2mat({ns_prop.EquivDiameter}'); %could also do
        %MajorAxisLength and MinorAxisLength, but length > width would mess
        %this up, so lets just assume length = 10 cm (as bins are 10cm
        %radius)
        list_width = (list_area)/10; % width in cm, assuming length is 10 cm
        list_width_max = max(list_width); % maximum width in radial segement (this is the "trunk" channel)
        chan_area18 = [chan_area18;list_area]; 
        chan_width18 = [chan_width18;list_width];
        chan_width_max18 = [chan_width_max18;list_width_max];
        n_chan = n;
        n_chan_i = [n_chan_i;n_chan];
    end
    area18(:,k) = mean(chan_area18, 'omitnan'); %cm^2
    area_std18(:,k) = std(chan_area18, 'omitnan'); %cm^2
    depth18(:,k) = mean(chan_depth18, 'omitnan'); %mm
    depth_std18(:,k) = std(chan_depth18, 'omitnan'); %mm
    width18(:,k) = mean(chan_width18, 'omitnan'); %cm
    width_std18(:,k) = std(chan_width18, 'omitnan'); %cm
    maxwidth18(:,k) = mean(chan_width_max18, 'omitnan'); %cm
    maxwidth_std18(:,k) = std(chan_width_max18, 'omitnan'); %cm
    n_chan18(:,k) = mean(n_chan_i, 'omitnan'); %number of channels
    n_chan_std18(:,k) = std(n_chan_i, 'omitnan'); %number of channels
end 

% treatment
% replace channel maps for treatment with NaN if time step does not have one
for i = 1:size(CM_19,3) 
    if sum(sum(CM_19(:,:,i), 'omitnan'), 'omitnan') == 0;
        CM_19(:,:,i) = NaN;
    end
end

% initialize empty matrices
area19 = []; % channel area
area_std19 =[]; % standard deviation of channel area
depth19 = []; % channel depth
depth_std19 = []; % standard deviation of channel depth 
width19 = []; % channel width
width_std19 = []; % standard deviation of channel width
maxwidth19 = []; % maximum channel width
maxwidth_std19 = []; % standard deviation of maximum channel width
n_chan19 = []; % number of channels
n_chan_std19 =[]; % standard deviation of the number of channels
dist = 0:100:3100; % distance we will calculate for (from 0 to 3100 mm, by 100 mm or 0.1 m).. we could change this if we want finer grained information

% run loop
for k = 1:(length(dist)-1)%loop to run through different radial distances from the end of the apex.
    %section to find x,y nodes for radial transect and generate matrix of
    %cross section topo and channel (yes/no) data
    idx = dd19(dd19 >= dist(k) & dd19 < dist(k+1));
    radial_dd = dd19 >= dist(k) & dd19 < dist(k+1);
    %how many channel pixels in the radius?
    chan_area19 = [];
    chan_depth19 = [];
    chan_width19 = [];
    chan_width_max19 = [];
    n_chan_i = [];
    for i = 1:size(CM_19,3) % loop through each timestep
        z = ZD_19(:,:,i); % elevation data for channel depth
        is_shot = CM_19(:,:,i).*radial_dd; % in channel or no?
        z = z.*is_shot; % elevation in the channels
        depth = [];
        if sum(sum(is_shot, 'omitnan')) == 0 % do not calculate channel properties for timesteps with no channel map
            n = NaN;
            ns_prop = NaN;
            list_area = NaN;
            list_width = NaN;
            list_width_max = NaN;
            depth = NaN;
        else %calculate channel properties for all other timesteps
            [label,n] = bwlabel(is_shot);
            ns_prop = regionprops(label, 'Area');
            list_area = (cell2mat({ns_prop.Area}'))*0.25; %cm^2
            %list_diam = cell2mat({ns_prop.EquivDiameter}'); %could also do
            %MajorAxisLength and MinorAxisLength, but length > width would mess
            %this up, so lets just assume length = 10 cm (as bins are 10cm
            %radius)
            % we will calculate depth for each channel here
            for j = 1:n % loop through all channels
                tmp = label;
                tmp(tmp ~= j) = NaN; % remove all data that is not in segement n
                tmp(tmp > 0) = 1; % turn segement data to 1
                tmp(tmp < 1) = NaN; % make everything else NaN
                z_tmp = z.*tmp; % segment elevations to get depth
                zs = z_tmp(~isnan(z_tmp)); % remove data not in the channel
                tmp_depth = (prctile(zs, 95) - prctile(zs, 5)); % channel depth in mm; we will take the 5% value as our channel bottom and the 95% value as our levee elevation, subtract the two.         
                depth = [depth, tmp_depth];
            end
            depth = max(depth); % we don't want to average in the non-trunk channels, as they are often a lot shallower
            list_width = (list_area)/10; %cm
            list_width_max = max(list_width);
        end
        chan_area19 = [chan_area19;list_area];
        chan_depth19 = [chan_depth19;depth];
        chan_width19 = [chan_width19;list_width];
        chan_width_max19 = [chan_width_max19;list_width_max];
        n_chan = n;
        n_chan_i = [n_chan_i;n_chan];
    end
    area19(:,k) = mean(chan_area19, 'omitnan'); % cm^2
    area_std19(:,k) = std(chan_area19, 'omitnan'); % cm^2
    depth19(:,k) = mean(chan_depth19, 'omitnan'); % mm
    depth_std19(:,k) = std(chan_depth19, 'omitnan'); % mm
    width19(:,k) = mean(chan_width19, 'omitnan'); % cm
    width_std19(:,k) = std(chan_width19, 'omitnan'); % cm
    maxwidth19(:,k) = mean(chan_width_max19, 'omitnan'); % cm
    maxwidth_std19(:,k) = std(chan_width_max19, 'omitnan'); % cm
    n_chan19(:,k) = mean(n_chan_i, 'omitnan'); % number of channels
    n_chan_std19(:,k) = std(n_chan_i, 'omitnan'); % number of channels
end 

%% Plot the data: Figure 3c
% create arrays
chan_array18 = [0:0.1:3; maxwidth18; maxwidth_std18];
cols = any(isnan(chan_array18),1);
chan_array18(:,cols) = [];

chan_array19 = [0:0.1:3; maxwidth19; maxwidth_std19];
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
fig = figure()
yyaxis left
patch([x18 fliplr(x18)], [curve1_18 fliplr(curve2_18)], 'b')
hold on
patch([x19 fliplr(x19)], [curve1_19 fliplr(curve2_19)], 'g')
alpha(0.15)
plot(x18, y18, 'b-', 'LineWidth', 2)
plot(x19, y19, 'g-', 'LineWidth', 2)
ylabel('max channel width (cm)')
ylim([-0.4 18])
yyaxis right
plot(0:0.1:3, n_chan18, 'bx')
plot(0:0.1:3, n_chan19, 'gx')
ylim([-0.1 4])
ylabel('number of channels')
xlabel('distance from apex (m)') 
legend('control mean', 'treatment mean', 'control stdev', 'treatment stdev', 'number of channels control', 'number of channels treatment')
ax = gca;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'k';
set(gcf, 'PaperUnits', 'inches');
y_width=6;x_width=8
set(gcf, 'PaperPosition', [0 0 x_width y_width]);
saveas(fig, '../figures/esurf_Figure3c.pdf')

%% Now we need to calculate far-field and channel aggradation rates
% This step will allow us to calculate channel in-filling timescale
% since the control experiment has LiDAR every hour but the treatment is
% every other hour, we will change the control to match treat to avoid any
% potential saddler effects (specifically important for calcuating
% aggradation rates
ZD_18 = ZD_18(:,:,2:2:560);
CM_18 = CM_18(:,:,2:2:560);
%Reload the channel maps, since we will use time steps with no channel
% map
cd '../data'
load('CM_19.mat')
CM_19 = CM_19(:,:,2:2:end); % remove every other map, so we can start with hour 1
cd '../code'

%% Process the elevation data to remove extranneous values

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
chanMaps_buffer_18 = [];
for j = 1:nt_18 %go by 2 so no saddler effects here
    buffer = imdilate(CM_18(:,:,j), se); %buffer channel map by 10 pixels
    chanMaps_buffer_18(:,:,j) = buffer;
end 
% treatment
chanMaps_buffer_19 = [];
for j = 1:nt_19
    buffer = imdilate(CM_19(:,:,j), se);
    chanMaps_buffer_19(:,:,j) = buffer;
end 

%% Clean treatment data
% Replace timesteps with no channel maps with the channel map from the
% next time step for the treatment experiment
for i = (nt_19-1):-1:1; %I know 280 has a channel map, so I can start at 279 and replace with channel map that comes after, this will work for ones that have multiple no maps in a row
    if sum(sum(chanMaps_buffer_19(:,:,i), 'omitnan'), 'omitnan') == 0;
        chanMaps_buffer_19(:,:,i) = chanMaps_buffer_19(:,:,i+1);
    end
end

%% Calculate channel and far-field aggradation
% Far field maps (remove channel from topo): control
FF_18_buffer = ~chanMaps_buffer_18; % inverse of channel buffer
ZD_18_FF_buffer = FF_18_buffer.*z18; % elevation of far field
ZD_18_FF_buffer(ZD_18_FF_buffer == 0) = NaN; % turn channel area into NaN
dZD_18_FF_buffer = diff(ZD_18_FF_buffer,1,3); % difference the maps along time dimension (3)

% Far field maps (remove channel from topo): treatment
FF_19_buffer = ~chanMaps_buffer_19; % inverse of channel buffer
ZD_19_FF_buffer = FF_19_buffer.*z19; % elevation of far field
ZD_19_FF_buffer(ZD_19_FF_buffer == 0) = NaN; % turn channel area into NaN
dZD_19_FF_buffer = diff(ZD_19_FF_buffer,1,3); % difference the maps along time dimension (3)

% Channel maps (remove far field from topo): control
ZD_18_chan_buffer = chanMaps_buffer_18.*z18; % elevation of channels
ZD_18_chan_buffer(ZD_18_chan_buffer == 0) = NaN; % remove far field area 
dZD_18_chan_buffer = diff(ZD_18_chan_buffer,1,3); % difference the maps along time dimension (3)

% Channel maps (remove far field from topo): treatment
ZD_19_chan_buffer = chanMaps_buffer_19.*z19; % elevation of channels
ZD_19_chan_buffer(ZD_19_chan_buffer == 0) = NaN; % remove far field area
dZD_19_chan_buffer = diff(ZD_19_chan_buffer,1,3); % difference the maps along time dimension (3)

%% Radial distance loop to get aggradation rates as a function of distance from the entrance channel
% control
chanMaps_buffer_18 = chanMaps_buffer_18(:,:,1:279); % when you difference the maps, you will have one less time step 
FF_18_buffer = FF_18_buffer(:,:,1:279); % far field maps
% initialize empty matrices to fill 
chan_agg_mean_18 = []; % mean channel aggradation
chan_agg_std_18 = []; % channel aggradation standard deviation
ff_agg_mean_18 = []; % mean far field aggradation
ff_agg_std_18 = []; % far field aggradation standard deviation
d_agg_18 = []; % difference in mean aggradation between channel and ff
d_agg_std18 = []; % propagated error
infill_18 = []; % channel infilling time
infill_std_18 = []; % channel infilling time
% loop through radial distances 
for k = 1:(length(dist)-1) % loop to run through different radial distances from the end of the entrance channel.
    % section to find x,y nodes for radial transect and generate matrix of
    % cross section topo and channel (yes/no) data
    idx = dd18(dd18 >= dist(k) & dd18 < dist(k+1));
    radial_dd = dd18 >= dist(k) & dd18 < dist(k+1); % can look at this using imagesc(radial_dd) to visualize a 0.1 m radial transect
    % initialize empty matrices 
    c_agg18 = []; 
    f_agg18 = [];
    for i = 1:(size(CM_18,3)-1)
        % channel and ff data in radial transect
        dz_chan = dZD_18_chan_buffer(:,:,i).*radial_dd; % multiply by radial transect
        dz_chan(dz_chan == 0.) = NaN; %remove channel pixels not in radial_dd
        dz_ff = dZD_18_FF_buffer(:,:,i).*radial_dd; % multiply by radial transect
        dz_ff(dz_ff == 0.) = NaN; % remove ff pixels not in radial_dd
        % now we will calculate channel aggradation vs. far field
        % aggradation
        dz_tmp = dz_chan(~isnan(dz_chan(:)));
        c_agg18 = [c_agg18;dz_tmp]; % mm/2-hr
        dz_tmp = dz_ff(~isnan(dz_ff(:)));
        f_agg18 = [f_agg18;dz_tmp]; % mm/2-hr 
    end
    %difference in agg between channel and ff
    dagg = (mean(c_agg18) - mean(f_agg18))/2; % mm/hr
    dagg_std = sqrt((std(c_agg18)^2)+(std(f_agg18)^2)/2); % error propagation
    % channel infilling time
    infill = depth18(k)/(dagg); % mm/mm/hr = hrs
    infill_std = infill*(sqrt((dagg_std/dagg)^2+(depth_std18(k)/depth18(k))^2)); % error propagation
    
    % save data
    chan_agg_mean_18(:,k) = mean(c_agg18, 'omitnan'); % mm/2-hr
    chan_agg_std_18(:,k) = std(c_agg18, 'omitnan'); % mm/2-hr
    ff_agg_mean_18(:,k) = mean(f_agg18, 'omitnan'); % mm/2-hr
    ff_agg_std_18(:,k) = std(f_agg18, 'omitnan'); % mm/2-hr
    d_agg_18(:,k) = dagg; % mm/hr
    d_agg_std18(:,k) = dagg_std; % mm/hr
    infill_18(:,k) = infill; % hours
    infill_std_18(:,k) = infill_std; % hours 
end 

% treatment
chanMaps_buffer_19 = chanMaps_buffer_19(:,:,1:279); % when you difference the maps, you will have one less time step 
FF_19_buffer = FF_19_buffer(:,:,1:279); % far field maps
% initialize empty matrices to fill 
chan_agg_mean_19 = []; % mean channel aggradation
chan_agg_std_19 = []; % channel aggradation standard deviation
ff_agg_mean_19 = []; % mean far field aggradation
ff_agg_std_19 = []; % far field aggradation standard deviation
d_agg_19 = []; % difference in mean aggradation between channel and ff
d_agg_std19 = []; % propagated error
infill_19 = []; % channel infilling time
infill_std_19 = []; % channel infilling time
% loop through radial distances 
for k = 1:(length(dist)-1) % loop to run through different radial distances from the end of the entrance channel.
    % section to find x,y nodes for radial transect and generate matrix of
    % cross section topo and channel (yes/no) data
    idx = dd19(dd19 >= dist(k) & dd19 < dist(k+1));
    radial_dd = dd19 >= dist(k) & dd19 < dist(k+1); % can look at this using imagesc(radial_dd) to visualize a 0.1 m radial transect
    % initialize empty matrices 
    c_agg19 = []; 
    f_agg19 = [];
    for i = 1:(size(CM_19,3)-1)
        % channel and ff data in radial transect
        dz_chan = dZD_19_chan_buffer(:,:,i).*radial_dd; % multiply by radial transect
        dz_chan(dz_chan == 0.) = NaN; %remove channel pixels not in radial_dd
        dz_ff = dZD_19_FF_buffer(:,:,i).*radial_dd; % multiply by radial transect
        dz_ff(dz_ff == 0.) = NaN; % remove ff pixels not in radial_dd
        % now we will calculate channel aggradation vs. far field
        % aggradation
        dz_tmp = dz_chan(~isnan(dz_chan(:)));
        c_agg19 = [c_agg19;dz_tmp]; % mm/2-hr
        dz_tmp = dz_ff(~isnan(dz_ff(:)));
        f_agg19 = [f_agg19;dz_tmp]; % mm/2-hr 
    end
    %difference in agg between channel and ff
    dagg = (mean(c_agg19) - mean(f_agg19))/2; % mm/hr
    dagg_std = sqrt((std(c_agg19)^2)+(std(f_agg19)^2)/2); % error propagation
    % channel infilling time
    infill = depth19(k)/(dagg); % mm/mm/hr = hrs
    infill_std = infill*(sqrt((dagg_std/dagg)^2+(depth_std19(k)/depth19(k))^2)); % error propagation
    
    % save data
    chan_agg_mean_19(:,k) = mean(c_agg19, 'omitnan'); % mm/2-hr
    chan_agg_std_19(:,k) = std(c_agg19, 'omitnan'); % mm/2-hr
    ff_agg_mean_19(:,k) = mean(f_agg19, 'omitnan'); % mm/2-hr
    ff_agg_std_19(:,k) = std(f_agg19, 'omitnan'); % mm/2-hr
    d_agg_19(:,k) = dagg; % mm/hr
    d_agg_std19(:,k) = dagg_std; % mm/hr
    infill_19(:,k) = infill; % hours
    infill_std_19(:,k) = infill_std; % hours 
end 

%% Plot the data: Figure 5
% error bars for plotting
% channel aggradation
chan_array18 = [0:0.1:3; chan_agg_mean_18; chan_agg_std_18];
cols = any(isnan(chan_array18),1);
chan_array18(:,cols) = [];

chan_array19 = [0:0.1:3; chan_agg_mean_19; chan_agg_std_19];
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
ff_array18 = [0:0.1:3; ff_agg_mean_18; ff_agg_std_18];
cols = any(isnan(ff_array18),1);
ff_array18(:,cols) = [];

ff_array19 = [0:0.1:3; ff_agg_mean_19; ff_agg_std_19];
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
dd_array18 = [0:0.1:3; d_agg_18; d_agg_std18];
cols = any(isnan(dd_array18),1);
dd_array18(:,cols) = [];

dd_array19 = [0:0.1:3; d_agg_19; d_agg_std19];
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
plot(x18, y18, 'b-', 'LineWidth', 2)
hold on
plot(x19, y19, 'g-', 'LineWidth', 2)
plot(ff_x18, ff_y18, 'b--', 'LineWidth', 2)
plot(ff_x19, ff_y19, 'g--', 'LineWidth', 2)

% Figure 5b,c
dist = 0:0.1:3;
fig2 = figure()
subplot(1,2,1)
patch([dd_x18 fliplr(dd_x18)], [dd_curve1_18 fliplr(dd_curve2_18)], 'b--')
hold on
patch([dd_x19 fliplr(dd_x19)], [dd_curve1_19 fliplr(dd_curve2_19)], 'g--')
alpha(0.15)
plot(dd_x18, dd_y18, 'b-', 'LineWidth', 2)
plot(dd_x19, dd_y19, 'g-', 'LineWidth', 2)
ylabel('channel in-filling rate (mm/hr)')
xlabel('distance from apex (m)') 
legend('control mean', 'treatment mean', 'control stdev', 'treatment stdev') 
set(gca,'XMinorTick','on','YMinorTick','on')
subplot(1,2,2)
plot(dist, infill_18, 'b-', 'LineWidth', 2)
hold on
plot(dist, infill_19, 'g-', 'LineWidth', 2)
ylim([0 300])
ylabel('channel in-filling time (hrs)')
xlabel('distance from apex (m)')
legend('control', 'treatment')
set(gca,'XMinorTick','on','YMinorTick','on')
y_width=7.25;x_width=9.125
set(gcf, 'PaperPosition', [0 0 x_width y_width]);
saveas(fig2, '../figures/esurf_Figure5bc.pdf')