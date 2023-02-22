clear all; close all
%% The purpose of this script is to calculate channel properties in the control and treatment experiments. 
% These data are plotted in Figure 4 in Sanks et al.
% (2023) submitted to Earth Surface Dynamics

% We will also calculate channel and far-field aggradation rates in the basin using a channel buffer as a function of radial distance from the apex.
% Finally, we will use the radial channel depths and radial channel and
% far-field aggradation rates to calculate a channel in-filling time
% (hours).

% Note that we will use every other hour of data for the control experiment
% in order to avoid any potential Saddler effects in our aggradation
% calculations, since the treatment experiment only contains LiDAR scans
% every other hour.

%% Load data
cd '..\data'
load('ZD_18.mat') %control topography (mm)
load('ZD_19.mat') %treatment topography (mm)
load('CM_18.mat') %control channel maps (binary, 1s = channels)
load('CM_19.mat') %treatment channel maps (binary, 1s = channels)
cd '..\code'

%% Define and set parameters
% control
% parameters
nx_18 = size(ZD_18, 1); %number of x locations on map
ny_18 = size(ZD_18,2); % number of y locations on map
nt_18 = size(ZD_18,3); % number of time steps in data set
dt_18 = 2; % delta t of time steps (hr); one hour initially but we changed data to be every other hour
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

%% Now we need to calculate far-field and channel aggradation rates
% This step will allow us to calculate channel in-filling timescale
% since the control experiment has LiDAR every hour but the treatment is
% every other hour, we will change the control to match treat to avoid any
% potential saddler effects (specifically important for calcuating
% aggradation rates
ZD_18 = ZD_18(:,:,2:2:560);
CM_18 = CM_18(:,:,2:2:560);

%% Process the elevation data to remove extranneous values
% control
% remove data below -9 mm rsl
nt_18 = size(CM_18, 3);
z18 = []; %go by 2 here so no saddler effects when comparing 18 and 19
for i = 1:nt_18
    z = ZD_18(:,:,i);
    z_rsl = z - (baselevel_rr*i*dt_18+ocean_zero);
    z_rsl(z_rsl<0) = NaN;
    z18(:,:,i) = z_rsl + (baselevel_rr*i*dt_18+ocean_zero);
end 
% treatment
%remove data below -9 mm rsl
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

%% Lets calculate the distance from apex to each pixel in our matrices
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

%% Radial distance loop to get aggradation rates as a function of distance from the apex
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
dist = 0:100:3100; % distance we will calculate for (from 0 to 3100 mm, by 100 mm or 0.1 m).. we could change this if we want finer grained information
% loop through radial distances 
for k = 1:(length(dist)-1) % loop to run through different radial distances from the end of the apex.
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
for k = 1:(length(dist)-1) % loop to run through different radial distances from the end of the apex.
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

%% Now lets plot the data
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

dist = 0:0.1:3;
fig = figure()
subplot(1,2,1)
% channel agg
% patch([x18 fliplr(x18)], [curve1_18 fliplr(curve2_18)], 'b-')
% hold on
% patch([x19 fliplr(x19)], [curve1_19 fliplr(curve2_19)], 'g-')
% alpha(0.15)
% plot(x18, y18, 'b-', 'LineWidth', 2)
% plot(x19, y19, 'g-', 'LineWidth', 2)
% dd agg
patch([dd_x18 fliplr(dd_x18)], [dd_curve1_18 fliplr(dd_curve2_18)], 'b--')
hold on
patch([dd_x19 fliplr(dd_x19)], [dd_curve1_19 fliplr(dd_curve2_19)], 'g--')
alpha(0.15)
plot(dd_x18, dd_y18, 'b--', 'LineWidth', 2)
plot(dd_x19, dd_y19, 'g--', 'LineWidth', 2)
%ylim([-2 5])
ylabel('aggradation rate (mm/hr)')
xlabel('distance from apex (m)') 
legend('control aggradation difference', 'treatment aggradation difference', 'control stdev', 'treatment stdev') %, 'control far-field', 'treatment far-field', 'control stdev', 'treatment stdev')
subplot(1,2,2)
plot(dist, infill_18, 'b*')
hold on
plot(dist, infill_19, 'g*')
ylim([0 300])
ylabel('channel in-filling time (hrs)')
xlabel('distance from apex (m)')
legend('control', 'treatment')
set(gcf, 'PaperUnits', 'inches');
y_width=7.25;x_width=9.125
set(gcf, 'PaperPosition', [0 0 x_width y_width]);

% dist = 0:0.1:3; % for plotting in m
% figure()
% yyaxis left
% errorbar(chan_agg_mean_18, chan_agg_std_18, 'b-')
% hold on
% errorbar(chan_agg_mean_19, chan_agg_std_19, 'color', [0 0.5 0],'linestyle','-')
% errorbar(ff_agg_mean_18, ff_agg_std_18, 'color', [0.5843 0.8157 0.9882], 'linestyle', '--')
% errorbar(ff_agg_mean_18, ff_agg_std_18, 'g--')
% yyaxis right
% plot(dist, infill_18, 'b*')
% hold on
% plot(dist, infill_19, 'g*')
