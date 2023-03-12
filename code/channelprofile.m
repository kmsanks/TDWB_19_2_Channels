clear all; close all
%% The purpose of this script is to create a channel profile from the treatment and control experiments
% These data are plotted in Figure 1c, 1c
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
load('ZD_18.mat') % control topography (mm)
load('ZD_19.mat') % treatment topography from Sam (mm)
load('terr_area18.mat') % binary map of area > -9 mm rsl
load('terr_area19.mat') % binary map of area > -9 mm rsl
cd '../code'

%% Crop data to every two hours for both experiments to match ZD_19
% This step will allow us to calculate channel in-filling timescale
% since the control experiment has LiDAR every hour but the treatment is
% every other hour, we will change the control to match treat to avoid any
% potential saddler effects (specifically important for calcuating
% aggradation rates
ZD_18 = ZD_18(:,:,1:2:560); % hour 1, 3, etc.
ZD_19 = ZD_19(:,:,2:end); % remove first time step because that is = 0 and we want to start at t = 1

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
baselevel_rr = 0.5; % base level rise rate (mm/2-hr)
ocean_zero = 25.25; % ocean elevation at hour 1 of each experiment (mm), hour 0 = 25 mm

% get boundary of basin, so we don't have issues later
% control
basin18 = ZD_18(:,:,1);
basin18(basin18 == 0) = NaN; % everything outside basin is NaN
basin18(~isnan(basin18)) = 1;% everything inside basin is 1
% treatment
basin19 = ZD_19(:,:,1);
basin19(basin19 == 0) = NaN; % well shadow is NaN
basin19(~isnan(basin19)) = 1;% everything inside basin is 1 except for well

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

%% Elevation relative to sea level

% control
z18 = [];
for i= 1:nt_18
    %What is sea level at time i
    sl = ((i-1)*baselevel_rr)+ocean_zero; %first wet scan starts at hour 1
    elevationmask = ZD_18(:,:,i); 
    elevationmask(elevationmask == 0) = NaN;
    elevationmask_rslr = elevationmask - sl;
    z18(:,:,i) = elevationmask_rslr;
end

% treatment
z19 = [];
for i= 1:nt_19
    %What is sea level at time i
    sl = ((i-1)*baselevel_rr)+ocean_zero; %first scan is t = 2, which is hour 1
    elevationmask = ZD_19(:,:,i); 
    elevationmask(elevationmask == 0) = NaN;
    elevationmask_rslr = elevationmask - sl;
    z19(:,:,i) = elevationmask_rslr;
end

% section to find x,y nodes for radial transect and generate matrix of
% cross section topo and channel (yes/no) data
rad_dist = 0:5:3100; % each radial transect is 50 mm or 5 cm
k = 200; % radial transect 100*5 = 500 mm or 1m from entrance
idx = dd18(dd18 >= rad_dist(k) & dd18 < rad_dist(k+1));
rads = dd18 >= rad_dist(k) & dd18 < rad_dist(k+1); % can look at this using imagesc(radial_dd) to visualize a 0.1 m radial transect% initialize empty matrices 
rads = rads.*basin18; %only radial on delta
rads = double(rads); % turn into double
rads(rads==0)=NaN; % so we can turn 0s to ones
% find x and y coords
[y18, x18] = find(~isnan(rads));

idx = dd19(dd19 >= rad_dist(k) & dd19 < rad_dist(k+1));
rads = dd19 >= rad_dist(k) & dd19 < rad_dist(k+1); % can look at this using imagesc(radial_dd) to visualize a 0.1 m radial transect% initialize empty matrices 
rads = rads.*basin19; %only radial on delta
rads = double(rads); % turn into double
rads(rads==0)=NaN; % so we can turn 0s to ones
% find x and y coords
[y19, x19] = find(~isnan(rads));

% which time do we want to plot
t18 = 40;
t19 = 110;
% Compute elevation values along transect using improfile
elevation_values18 = improfile(z18(:,:,t18), x18, y18);
radial_dist18 = 0:0.005:(length(elevation_values18)-1)*0.005; %distance along radial transect

elevation_values19 = improfile(z19(:,:,t19), x19, y19);
radial_dist19 = 0:0.005:(length(elevation_values19)-1)*0.005; %distance along radial transect

cmap = colormap("parula");

imAlpha=ones(size(z18(:,:,t18)));
imAlpha(isnan(z18(:,:,t18)))=0;
imagesc(z18(:,:,t18),'AlphaData',imAlpha);
set(gca,'color',0*[1 1 1]);

% Add topography data as a contour plot
fig = figure;
imagesc(z18(:,:,t18), 'AlphaData', ~isnan(z18(:,:,t18))) % Plot the topography data as a filled contour plot with 50 levels
axis equal
xlim([0 500])
ylim([109 500])
hold on
plot(x18, y18, 'k-', 'LineWidth', 2) % Plot the radial line in black with a dashed line style
colormap(cmap)
colorbar
clim([-9 40]) % Add a colorbar to show the elevation scale
legend('radial transect (1e)')
set(gca, 'XMinorTick', 'On', 'YMinorTick', 'On')
set(gcf, 'PaperUnits', 'inches');
y_width=7.25 ;x_width=9.125;
set(gcf, 'PaperPosition', [0 0 x_width y_width]);
saveas(gcf,'../figures/esurf_Figure1c.pdf')

fig = figure;
imagesc(z19(:,:,t19), 'AlphaData', ~isnan(z19(:,:,t19))) % Plot the topography data as a filled contour plot with 50 levels
axis equal
xlim([130 630])
ylim([214 605])
hold on
plot(x19, y19, 'k-', 'LineWidth', 2) % Plot the radial line in black with a dashed line style
colormap(cmap)
colorbar % Add a colorbar to show the elevation scale
clim([-19 20])
legend('radial transect (1f)')
set(gca, 'XMinorTick', 'On', 'YMinorTick', 'On')
set(gcf, 'PaperUnits', 'inches');
y_width=7.25 ;x_width=9.125;
set(gcf, 'PaperPosition', [0 0 x_width y_width]);
saveas(gcf,'../figures/esurf_Figure1d.pdf')

fig = figure;
plot(radial_dist18, elevation_values18, 'b-', 'LineWidth', 1)
xlabel('distance along radial transect (m)')
ylabel('elevation relative to sea level (mm)')
ylim([2 33])
xlim([0 2.2])
pbaspect([4 0.5 1])
set(gca, 'XMinorTick', 'On', 'YMinorTick', 'On')
set(gcf, 'PaperUnits', 'inches');
y_width=7.25 ;x_width=9.125;
set(gcf, 'PaperPosition', [0 0 x_width y_width]);
saveas(gcf,'../figures/esurf_Figure1e.pdf')

fig = figure;
plot(radial_dist19, elevation_values19, 'g-', 'LineWidth', 1)
xlabel('distance along radial transect (m)')
ylabel('elevation relative to sea level (mm)')
ylim([-10 25])
xlim([0 2.2])
pbaspect([4 0.5 1])
set(gca, 'XMinorTick', 'On', 'YMinorTick', 'On')
set(gcf, 'PaperUnits', 'inches');
y_width=7.25 ;x_width=9.125;
set(gcf, 'PaperPosition', [0 0 x_width y_width]);
saveas(gcf,'../figures/esurf_Figure1f.pdf')