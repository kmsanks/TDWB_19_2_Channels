clear all; close all
%% The purpose of this script is to calculate channel and far-field aggradation and channel in-filling rate as a function of radial distance
% These data are plotted in Table 1 and Figure 5a in Sanks et al.
% (2023) submitted to Earth Surface Dynamics

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
ZD_18 = ZD_18(:,:,2:2:560); % go by 2 to avoid Saddler
CM_18 = CM_18(:,:,2:2:560); % go by 2 to avoid Saddler
nx_18 = size(ZD_18, 1); %number of x locations on map
ny_18 = size(ZD_18,2); % number of y locations on map
nt_18 = size(ZD_18,3); % number of time steps in data set
dt_18 = 2; % delta t of time steps (hr)
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

%% Remove data below -9 mm rsl
% control
% remove data below -9 mm rsl
z18 = []; %go by 2 here so no saddler effects when comparing 18 and 19
for i = 1:nt_18 % loop through time
    z = ZD_18(:,:,i); % elevations
    z_rsl = z - (baselevel_rr*(i+1)*dt_18+ocean_zero); % since we removed hour 1, i+1
    z_rsl(z_rsl<-9) = NaN; % data below -9 = NaN
    z18(:,:,i) = z_rsl + (baselevel_rr*(i+1)*dt_18+ocean_zero); % return back to not relative to sea level
end 

% treatment
%remove data below -9 mm rsl
z19 = [];
for i = 1:nt_19 % loop through time
    z = ZD_19(:,:,i); % elevations
    z_rsl = z - (baselevel_rr*i*dt_19+ocean_zero); %since we start on timestep 2, which is hour 1 use i
    z_rsl(z_rsl<-9) = NaN; % data below -9 = NaN
    z19(:,:,i) = z_rsl + (baselevel_rr*i*dt_19+ocean_zero); %return back to not relative to sea level
end 

%% Buffer the channels

% control
se = strel('square',10); %create a square of 10 pixels around each 1 point, this will buffer our channel maps
CM_buffer_18 = [];
for j = 1:nt_18 
    buffer = imdilate(CM_18(:,:,j), se); % buffer channel map by 10 pixels
    CM_buffer_18(:,:,j) = buffer;
end 

% treatment
CM_buffer_19 = [];
for j = 1:nt_19
    buffer = imdilate(CM_19(:,:,j), se); % buffer channel map by 10 pixels
    CM_buffer_19(:,:,j) = buffer;
end 

%% Replace timesteps with no channel maps with the channel map from the next time step for the treatment experiment
for i = (nt_19-1):-1:1; %I know 280 has a channel map, so I can start at 279 and replace with channel map that comes after, this will work for ones that have multiple no maps in a row
    if sum(sum(CM_buffer_19(:,:,i), 'omitnan'), 'omitnan') == 0
        CM_buffer_19(:,:,i) = CM_buffer_19(:,:,i+1);
    end
end

%% Calculate far field and channel aggradation rates

% far field
% control
% far field maps (remove channel from topo)
FF_18_buffer = ~CM_buffer_18;
ZD_18_FF_buffer = FF_18_buffer.*z18;
ZD_18_FF_buffer(ZD_18_FF_buffer == 0) = NaN;
% elevation difference
dZD_18_FF_buffer = diff(ZD_18_FF_buffer,1,3);

% treatment
% far field maps (remove channel from topo)
FF_19_buffer = ~CM_buffer_19;
ZD_19_FF_buffer = FF_19_buffer.*z19;
ZD_19_FF_buffer(ZD_19_FF_buffer == 0) = NaN;
%elevation difference
dZD_19_FF_buffer = diff(ZD_19_FF_buffer,1,3);

% channels
% control
% remove far field agg from topo (we only want the channels)
ZD_18_chan_buffer = CM_buffer_18.*z18;
ZD_18_chan_buffer(ZD_18_chan_buffer == 0) = NaN;
% elevation difference
dZD_18_chan_buffer = diff(ZD_18_chan_buffer,1,3);

% treatment
% remove far field agg from topo (we only want the channels)
ZD_19_chan_buffer = CM_buffer_19.*z19;
ZD_19_chan_buffer(ZD_19_chan_buffer == 0) = NaN;
% elevation difference
dZD_19_chan_buffer = diff(ZD_19_chan_buffer,1,3);

%% Calculate aggradation statistics for Table 1
mean_chan_agg18 = mean(dZD_18_chan_buffer(:), 'omitnan')/2; %mm/hr
stdev_chan_agg18 = std(dZD_18_chan_buffer(:), 'omitnan')/2; %mm/hr

mean_chan_agg19 = mean(dZD_19_chan_buffer(:), 'omitnan')/2; %mm/hr
stdev_chan_agg19 = std(dZD_19_chan_buffer(:), 'omitnan')/2; %mm/hr

mean_ff_agg18 = mean(dZD_18_FF_buffer(:), 'omitnan')/2; %mm/hr
stdev_ff_agg18 = std(dZD_18_FF_buffer(:), 'omitnan')/2; %mm/hr

mean_ff_agg19 = mean(dZD_19_FF_buffer(:), 'omitnan')/2; %mm/hr
stdev_ff_agg19 = std(dZD_19_FF_buffer(:), 'omitnan')/2; %mm/hr

clear('CM_18', 'CM_19', 'z18', 'z19', 'z_rsl');

%% Calculate aggradation rates and channel in-filling rate as a function of radial distance
% control
CM_buffer_18 = CM_buffer_18(:,:,1:279);
FF_18_buffer = FF_18_buffer(:,:,1:279);
chan_agg_mean_18 = [];
chan_agg_std_18 = [];
ff_agg_mean_18 = [];
ff_agg_std_18 = [];
ns_chan_18 = [];
ns_ff_18 = [];
for i = 1:650 %loop to run through different radial distances from the end of the entrance channel.
    %section to find x,y nodes for radial transect and generate matrix of
    %cross section topo and channel (yes/no) data
    th = 0:1/i:2*pi;
    xunit = i * cos(th) + xentrance_18;
    yunit = i * sin(th) + yentrance_18;
    xunit = round(xunit);
    yunit = round(yunit);
    xunit18 = [];
    yunit18 = [];
    chandz_18_list = [];
    ffz_18_list = [];
    for j = 1:max(size(xunit))
        xc = xunit(j);
        yc = yunit(j);
        if xc >= 1;
            if yc >= 1;
                if xc <= nx_18
                    if yc <= ny_18
                        xunit18 = [xunit18;xc];
                        yunit18 = [yunit18;yc];
                    end
                end
            end
        end
    end
    xunit = xunit18;
    yunit = yunit18;
    XY18 = [yunit xunit];
    XY18 = sortrows(XY18);
    xunit = XY18(:,2);
    yunit = XY18(:,1);
    zs18 = [];
    is18 = [];
    for j = 1:max(size(xunit));
        xs_shot = dZD_18_chan_buffer(xunit(j),yunit(j),:);
        is_shot = CM_buffer_18(xunit(j),yunit(j),:);
        zs18 = [zs18;xs_shot];
        is18 = [is18;is_shot];
 
    end
    zs18 = squeeze(zs18);%matrix of topo along radial transect
    is18 = squeeze(is18);%matrix of channel (yes/no) along radial transect
    
    %find channel dz
    [chan_row, chan_col] = find(is18 == 1);
    chandz_18_list = zs18(sub2ind(size(zs18), chan_row, chan_col));
    ns_chan = length(chandz_18_list(~isnan(chandz_18_list))); %number of all non-nan values
    chan_agg_mean = mean(chandz_18_list, 'omitnan');
    chan_agg_std = std(chandz_18_list, 'omitnan');
    chan_agg_mean_18 = [chan_agg_mean_18; chan_agg_mean];
    chan_agg_std_18 = [chan_agg_std_18; chan_agg_std];
    ns_chan_18 = [ns_chan_18; ns_chan];
    
    zs18 = [];
    is18 = [];
    for j = 1:max(size(xunit));
        xs_shot = dZD_18_FF_buffer(xunit(j),yunit(j),:);
        is_shot = FF_18_buffer(xunit(j),yunit(j),:);
        zs18 = [zs18;xs_shot];
        is18 = [is18;is_shot];
    end
    zs18 = squeeze(zs18);%matrix of topo along radial transect
    is18 = squeeze(is18);%matrix of channel (yes/no) along radial transect
    
    %find channel dz
    [ff_row, ff_col] = find(is18 == 1);
    ffdz_18_list = zs18(sub2ind(size(zs18), ff_row, ff_col));
    ns_ff = length(ffdz_18_list(~isnan(ffdz_18_list))); %number of all non-nan values
    ff_agg_mean = mean(ffdz_18_list, 'omitnan');
    ff_agg_std = std(ffdz_18_list, 'omitnan');
    ff_agg_mean_18 = [ff_agg_mean_18; ff_agg_mean];
    ff_agg_std_18 = [ff_agg_std_18; ff_agg_std];
    ns_ff_18 = [ns_ff_18; ns_ff];
end

% treatment
CM_buffer_19 = CM_buffer_19(:,:,1:279);
FF_19_buffer = FF_19_buffer(:,:,1:279);
chan_agg_mean_19 = [];
chan_agg_std_19 = [];
ff_agg_mean_19 = [];
ff_agg_std_19 = [];
ns_chan_19 = [];
ns_ff_19 = [];
%%look to see if agg rate is divided by dt
for i = 1:650;%loop to run through different radial distances from the end of the entrance channel.
    %section to find x,y nodes for radial transect and generate matrix of
    %cross section topo and channel (yes/no) data
    th = 0:1/i:2*pi;
    xunit = i * cos(th) + xentrance_19;
    yunit = i * sin(th) + yentrance_19;
    xunit = round(xunit);
    yunit = round(yunit);
    xunit19 = [];
    yunit19 = [];
    chandz_19_list = [];
    ffz_19_list = [];
    for j = 1:max(size(xunit))
        xc = xunit(j);
        yc = yunit(j);
        if xc >= 1;
            if yc >= 1;
                if xc <= nx_19
                    if yc <= ny_19
                        xunit19 = [xunit19;xc];
                        yunit19 = [yunit19;yc];
                    end
                end
            end
        end
    end
    xunit = xunit19;
    yunit = yunit19;
    XY19 = [yunit xunit];
    XY19 = sortrows(XY19);
    xunit = XY19(:,2);
    yunit = XY19(:,1);
    zs19 = [];
    is19 = [];
    for j = 1:max(size(xunit));
        xs_shot = dZD_19_chan_buffer(xunit(j),yunit(j),:);
        is_shot = CM_buffer_19(xunit(j),yunit(j),:);
        zs19 = [zs19;xs_shot];
        is19 = [is19;is_shot];
    end
    zs19 = squeeze(zs19);%matrix of topo along radial transect
    is19 = squeeze(is19);%matrix of channel (yes/no) along radial transect
    
    %find channel dz
    [chan_row, chan_col] = find(is19 == 1);
    if length(chan_row) == 0
        chandz_19_list = NaN;
    else 
        chandz_19_list = zs19(sub2ind(size(zs19), chan_row, chan_col));
    end
    ns_chan = length(chandz_19_list(~isnan(chandz_19_list))); %number of all non-nan values
    chan_agg_mean = mean(chandz_19_list, 'omitnan');
    chan_agg_std = std(chandz_19_list, 'omitnan');
    chan_agg_mean_19 = [chan_agg_mean_19; chan_agg_mean];
    chan_agg_std_19 = [chan_agg_std_19; chan_agg_std];
    ns_chan_19 = [ns_chan_19; ns_chan];
    
    zs19 = [];
    is19 = [];
    for j = 1:max(size(xunit));
        xs_shot = dZD_19_FF_buffer(xunit(j),yunit(j),:);
        is_shot = FF_19_buffer(xunit(j),yunit(j),:);
        zs19 = [zs19;xs_shot];
        is19 = [is19;is_shot];
    end
    zs19 = squeeze(zs19);%matrix of topo along radial transect
    is19 = squeeze(is19);%matrix of channel (yes/no) along radial transect
    
    %find channel dz
    [ff_row, ff_col] = find(is19 == 1);
    ffdz_19_list = zs19(sub2ind(size(zs19), ff_row, ff_col));
    ns_ff = length(ffdz_19_list(~isnan(ffdz_19_list))); %number of all non-nan values
    ff_agg_mean = mean(ffdz_19_list, 'omitnan');
    ff_agg_std = std(ffdz_19_list, 'omitnan');
    ff_agg_mean_19 = [ff_agg_mean_19; ff_agg_mean];
    ff_agg_std_19 = [ff_agg_std_19; ff_agg_std];
    ns_ff_19 = [ns_ff_19; ns_ff];
end

%% Calculate channel in-filling rate for Table 1
chan_mean_18 = mean(chan_agg_mean_18, 'omitnan')/2;
ff_mean_18 = mean(ff_agg_mean_18, 'omitnan')/2;
chan_std_18 = mean(chan_agg_std_18, 'omitnan')/2;
ff_std_18 = mean(ff_agg_std_18, 'omitnan')/2;

chan_mean_19 = mean(chan_agg_mean_19, 'omitnan')/2;
ff_mean_19 = mean(ff_agg_mean_19, 'omitnan')/2;
chan_std_19 = mean(chan_agg_std_19, 'omitnan')/2;
ff_std_19 = mean(ff_agg_std_19, 'omitnan')/2;

% channel in-filling rate
% ********** NEED TO FIX WITH CHANNEL DEPTH ISSUE ***********************
%control
c_if18 = 14.6/(chan_mean_18 - ff_mean_18); %mm/(mm/hr) = hrs %channel depth/(in-channel agg - ff_agg)
c_if_add18 = sqrt((chan_std_18)^2+(ff_std_18)^2);
c_if_se18 = c_if18*sqrt((c_if_add18/(chan_mean_18 - ff_mean_18))^2 + (7.82/14.6)^2);

% treatment
c_if19 = 13.4/(chan_mean_19 - ff_mean_19); %mm/(mm/hr) = hrs
c_if_add19 = sqrt((chan_std_19)^2+(ff_std_19)^2);
c_if_se19 = c_if19*sqrt((c_if_add19/(chan_mean_19 - ff_mean_19))^2 + (6.03/13.4)^2);

%% Plot the data for Figure 5a
dist_18 = 1:650; % radial distances 
dist_19 = 1:650; % radial distances
fig = figure();
%yyaxis left
plot((dist_18*5)/1000, chan_agg_mean_18', 'b-', 'linewidth',2);
hold on
plot((dist_19*5)/1000, chan_agg_mean_19', 'g-', 'linewidth',2);
plot((dist_18*5)/1000, ff_agg_mean_18', 'b--', 'linewidth',2);
plot((dist_19*5)/1000, ff_agg_mean_19', 'g--', 'linewidth',2);
yline(0.50,  'k:', 'linewidth', 2);
ylabel('mean aggradation rate (mm/2-hr)');
ylim([0 4]);
xlim([0 3.5]);
xlabel('radial distance from apex (m)');
legend('control in-channel', 'treatment in-channel',...
    'control far-field', 'treatment far-field',...
    'RSLRb')
set(gca,'XMinorTick','on','YMinorTick','on')
saveas(fig, '../figures/esurf_Figure5a.pdf');

% Look at if number of pixels impacts these aggradation rates
figure()
subplot(2,1,1)
plot((dist_18*5)/1000, chan_agg_mean_18', 'b-', 'linewidth',2);
hold on
plot((dist_19*5)/1000, chan_agg_mean_19', 'g-', 'linewidth',2);
plot((dist_18*5)/1000, ff_agg_mean_18', 'b--', 'linewidth',2);
plot((dist_19*5)/1000, ff_agg_mean_19', 'g--', 'linewidth',2);
yline(0.25,  'k:', 'linewidth', 2);
ylabel('mean aggradation rate (mm/hr)');
legend('control in-channel', 'treatment in-channel',...
    'control far-field', 'treatment far-field',...
    'RSLRb')
set(gca,'XMinorTick','on','YMinorTick','on')
subplot(2,1,2)
plot((dist_18*5)/1000, ns_chan_18', 'b-', 'linewidth', 1);
hold on
plot((dist_19*5)/1000, ns_chan_19', 'g-', 'linewidth', 1);
plot((dist_18*5)/1000, ns_ff_18', 'b--', 'linewidth', 1);
plot((dist_19*5)/1000, ns_ff_19', 'g--', 'linewidth', 1);
yline(100, 'k-');
ylim([0 5000])
legend('control n chan', 'treatment n chan',...
     'control n far-field', 'treatment n far-field', 'n = 100');
xlabel('radial distance from the entrance (m)');
ylabel('number of pixels');
