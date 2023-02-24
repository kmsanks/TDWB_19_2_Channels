clear all; close all
%% The purpose of this script is to calculate channel bed elevation as a function of radial distance from the apex
% These data are plotted in Figure 3a in Sanks et al.
% (2023) submitted to Earth Surface Dynamics
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

%% Extract channel topo
% Make elevation screen, so flow can be clipped to the areas above sea level

% control
tmp18 = [];
for i= 1:nt_18
    %What is sea level at time i
    sl = (i*baselevel_rr*dt_18)+25; %first wet scan starts 48 minutes into the hour so i = i is a good approximation
    elevationmask = ZD_18(:,:,i); 
    elevationmask(elevationmask == 0) = NaN;
    elevationmask_rslr = elevationmask - sl;
    %elevationmask_rslr(elevationmask_rslr < -30) = NaN; %
    tmp18(:,:,i) = elevationmask_rslr;
end

% treatment
tmp19 = [];
for i= 1:nt_19
    %What is sea level at time i
    sl = ((i-1)*baselevel_rr*dt_19)+25; %first scan is at hour 0, so we need i-1
    elevationmask = ZD_19(:,:,i); 
    elevationmask(elevationmask == 0) = NaN;
    elevationmask_rslr = elevationmask - sl;
    %elevationmask_rslr(elevationmask_rslr < -30) = NaN;
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

%% Loop to find channel bed elevation relative to sea level as a function of radial distance from apex

elev_18 = [];%empty array to be filled with elevation from radial transects every 5 mm from source 
elev_std_18 = [];%empty array to be filled with std elevation from radial transects every 5 mm from source 
elev_25prct_18 = [];%fill loop for 1.5 iqr
elev_75prct_18 = [];

%Loop through radial transects
for i = 1:650%loop to run through different radial distances from the end of the entrance channel.
    %section to find x,y nodes for radial transect and generate matrix of
    %cross section topo and channel (yes/no) data
    th = 0:1/i:2*pi;
    xunit = i * cos(th) + xentrance_18;
    yunit = i * sin(th) + yentrance_18;
    xunit = round(xunit);
    yunit = round(yunit);
    xunit18 = [];
    yunit18 = [];
    for j = 1:max(size(xunit))
        xc = xunit(j);
        yc = yunit(j);
        if xc >= 1
            if yc >= 1
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
    for j = 1:max(size(xunit))
        xs_shot = z18(xunit(j),yunit(j),:);
        zs18 = [zs18;xs_shot];
    end
    zs18 = squeeze(zs18);%matrix of topo along radial transect

    elev_18(:,i) = mean(zs18, 'omitnan');%empty array to be filled with elevation from radial transects
    elev_25prct_18(:,i) = prctile(zs18,25);
    elev_25prct_18(:,i) = prctile(zs18,25);
    elev_std_18(:,i) = std(zs18, 'omitnan');
    radial_dist_18(:,i) = i*.005; %m
end

elev_18_mean = mean(elev_18, 'omitnan');
elev_18_stdmean = std(elev_18, 'omitnan');

elev_19 = [];%empty array to be filled with elevation from radial transects every 5 mm from source 
elev_std_19 = [];%empty array to be filled with std elevation from radial transects every 5 mm from source 
elev_25prct_19 = [];
elev_75prct_19 = [];

%Loop through radial transects
for i = 1:666%loop to run through different radial distances from the end of the entrance channel.
    %section to find x,y nodes for radial transect and generate matrix of
    %cross section topo and channel (yes/no) data
    th = 0:1/i:2*pi;
    xunit = i * cos(th) + xentrance_19;
    yunit = i * sin(th) + yentrance_19;
    xunit = round(xunit);
    yunit = round(yunit);
    xunit19 = [];
    yunit19 = [];
    for j = 1:max(size(xunit))
        xc = xunit(j);
        yc = yunit(j);
        if xc >= 1
            if yc >= 1
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
    for j = 1:max(size(xunit))
        xs_shot = z19(xunit(j),yunit(j),:);
        zs19 = [zs19;xs_shot];
    end
    zs19 = squeeze(zs19);%matrix of topo along radial transect

    elev_19(:,i) = mean(zs19, 'omitnan');%empty array to be filled with elevation from radial transects
    elev_25prct_19(:,i) = prctile(zs19,25);
    elev_75prct_19(:,i) = prctile(zs19,75);
    elev_std_19(:,i) = std(zs19, 'omitnan');
    radial_dist_19(:,i) = i*.005; %m
end 

elev_19_mean = mean(elev_19, 'omitnan');
elev_19_stdmean = std(elev_19, 'omitnan');

%% Plot the data
ybars = [-9 5];

fig = figure()
plot(radial_dist_18/1000, elev_18_mean, 'b-', 'linewidth', 2)
hold on
plot(radial_dist_19/1000, elev_19_mean, 'g-', 'linewidth', 2)
yline(0, 'k-', 'linewidth', 2)
patch([min(xlim) max(xlim) max(xlim) min(xlim)], [ybars(1) ybars(1), ybars(2) ybars(2)], [0.8 0.8 0.8])
plot(radial_dist_18/1000, elev_18_mean, 'b-', 'linewidth', 2)
plot(radial_dist_19/1000, elev_19_mean, 'g-', 'linewidth', 2)
yline(0, 'k-', 'linewidth', 2)
%yline(-9, 'r-')
%yline(5, 'r-')
xlabel('radial distance from entrance (mm)')
ylabel('mean channel bed elevation relative to sea level (mm)')
legend('control', 'treatment', 'sea level', 'marsh window')
set(gcf, 'PaperUnits', 'inches');
y_width=7.25;x_width=9.125
set(gcf, 'PaperPosition', [0 0 x_width y_width]);
saveas(fig, 'radialelevation.pdf')

array18 = [radial_dist_18/1000; elev_18_mean; elev_18_stdmean];
cols = any(isnan(array18),1);
array18(:,cols) = [];

array19 = [radial_dist_19/1000; elev_19_mean; elev_19_stdmean];
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

fig2 = figure()
plot(x18, y18, 'b', 'LineWidth', 2)
hold on
plot(x19, y19, 'g', 'LineWidth', 2)
patch([x18 fliplr(x18)], [curve1_18 fliplr(curve2_18)], 'b')
patch([x19 fliplr(x19)], [curve1_19 fliplr(curve2_19)], 'g')
%yline(-9, 'k-')
%yline(5, 'k-')
patch([min(xlim) max(xlim) max(xlim) min(xlim)], [ybars(1) ybars(1), ybars(2) ybars(2)], 'k')
alpha(0.15)
yline(0, 'k-', 'linewidth', 2)
plot(x18, y18, 'b', 'LineWidth', 2)
plot(x19, y19, 'g', 'LineWidth', 2)
ylim([-40 50])
xlim([0 3])
ylabel('mean channel bed elevation relative to sea level (mm)')
xlabel('radial distance from channel entrance (m)') 
legend('control mean', 'treatment mean', 'control stdev', 'treatment stdev', 'marsh window', 'sea level')
%patch([min(xlim) max(xlim) max(xlim) min(xlim)], [ybars(1) ybars(1), ybars(2) ybars(2)], [0.8 0.8 0.8])
set(gcf, 'PaperUnits', 'inches');
y_width=7.25;x_width=9.125
set(gcf, 'PaperPosition', [0 0 x_width y_width]);
saveas(fig2, 'radial_channel_elevation.pdf')
