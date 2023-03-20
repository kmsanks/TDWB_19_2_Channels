clear all; close all;
%% The purpose of this script is to calculate lateral mobility statistics 
% This script creates Figure 6 and SI Figs. B1, B2, B6, and B7, as well as some results from Table 1 presented in Sanks et al.(2023) submitted to Earth Surface Dynamics

%% First, lets set our relative path
cd '../' %We should have opened this code from 'C:\Users\XXXXXXX\Documents\GitHub\TDWB_19_2_Channels\code'

%% Load data
cd '.\data'
load('ZD_18.mat');
load('CM_18.mat');
load('ZD_19.mat');
load('CM_19.mat');
cd '..\code'

%% Set parameters
%control
nx_18 = size(ZD_18,1); %number of x locations on map
ny_18 = size(ZD_18,2); %number of y locations on map
nt_18 = size(ZD_18,3); %number of time steps in data set
dt_18 = 1; %delta t of time steps (hr)
xentrance_18 = 109; %x grid node location of the entrance channel
yentrance_18 = 271; %y grid node location of the entrance channel

%treatment
nx_19 = size(ZD_19,1); %number of x locations on map
ny_19 = size(ZD_19,2); %number of y locations on map
nt_19 = size(ZD_19,3); %number of time steps in data set
dt_19 = 2; %delta t of time steps (hr)
xentrance_19 = 214; %x grid node location of the entrance channel (x is down dip)
yentrance_19 = 397; %y grid node location of the entrance channel (y is strike)

%both experiments
baselevel_rr = 0.25; %base level rise rate (mm/hr)
ocean_zero = 25; %ocean elevation at beginning of experiment (mm)
dx = 5; %mm grid cells

% get boundary of basin, so we don't have issues later
% control
basin18 = ZD_18(:,:,1);
basin18(basin18 == 0) = NaN; % everything outside basin is NaN
basin18(~isnan(basin18)) = 1;% everything inside basin is 1
% treatment
basin19 = ZD_19(:,:,1);
basin19(basin19 == 0) = NaN; % well shadow is NaN
basin19(~isnan(basin19)) = 1;% everything inside basin is 1 except for well

%% Clean treatment channel data
%replace empty channel maps in treatment with channel map during timestep
%after
for i = (nt_19-1):-1:1 %I know 560 has a channel map, so I can start at 279 and replace with channel map that comes after, this will work for ones that have multiple no maps in a row
    if sum(sum(CM_19(:,:,i), 'omitnan'), 'omitnan') == 0
        CM_19(:,:,i) = CM_19(:,:,i+1);
    end
end

% %% Find area above sea level for at least 50% of the experiment
% % control
% pland_18 = [];
% t_18 = [];
% for i = 1:nt_18
%     t_18 = [t_18;i*dt_18]; %time (hours)
%     z = ZD_18(:,:,i); %elevation (mm)
%     z = z - (baselevel_rr*i*dt_18+ocean_zero); %elevation relative to sea level
%     z(z > 0) = 1; %above sea level binary of 1s
%     z(z < 0) = 0; %below sea level binary of 0s
%     pland_18(:,:,i) = z; %save land mask
% end
% 
% % treatment
% znt_19 = size(ZD_19,3);
% zdt_19 = 2;
% pland_19 = [];
% t_19 = [];
% for i = 1:znt_19
%     t_19 = [t_19;(i-1)*zdt_19]; %time (hours)
%     z = ZD_19(:,:,i); %elevation (mm)
%     z = z - (baselevel_rr*(i-1)*zdt_19+ocean_zero); %elevation relative to sea level
%     z(z > 0) = 1; %above sea level binary of 1s
%     z(z < 0) = 0; %below sea level binary of 0s
%     pland_19(:,:,i) = z; %save land mask
% end 
% 
% %control
% pland_18_sum = sum(pland_18,3); %how many times each grid cell is above sea level
% pland_18_sum(pland_18_sum == 0) = NaN; %not on land make NaN
% pland_18_frac = pland_18_sum/560; %fraction of time pixel is on land
% 
% frac_18 = [];
% area_pland_18 = []; %area on land for various amounts of the experiment
% for i = round(linspace(0,1,101),2)
%     frac_18 = [frac_18;i];
%     area = pland_18_frac(pland_18_frac >= i); %& pland_18_frac < i+1
%     area(area > 0) = 1;
%     areasum = sum(area);
%     area_pland_18_i = (areasum)*2.5*10^-5; %area of 1 pixel
%     area_pland_18 = [area_pland_18; area_pland_18_i];     
% end 
% 
% %treatment
% pland_19_sum = sum(pland_19,3);
% pland_19_sum(pland_19_sum == 0) = NaN;
% pland_19_frac = pland_19_sum/281;
% 
% frac_19 = [];
% area_pland_19 = [];
% for i = round(linspace(0,1,101),2)
%     frac_19 = [frac_19;i];
%     area = pland_19_frac(pland_19_frac >= i); %& pland_18_frac < i+1
%     area(area > 0) = 1;
%     areasum = sum(area);
%     area_pland_19_i = (areasum)*2.5*10^-5; %area of 1 pixel
%     area_pland_19 = [area_pland_19; area_pland_19_i];     
% end 
% 
% %% Crop channel maps to 50% land area
% % crop channels to 50% land area
% pland_18_frac_nan = pland_18_frac;
% pland_18_frac_nan(pland_18_frac_nan < 0.5) = NaN;
% pland_18_frac_nan(pland_18_frac_nan > 0) = 1;
% CM_18_crop = [];
% for i = 1:((size(CM_18,3)))
%     CM_18_nan = CM_18(:,:,i).*pland_18_frac_nan; % crop channel maps to 50% land area
%     CM_18_crop(:,:,i) = CM_18_nan; % save new channel maps
% end
% 
% % treatment
% pland_19_frac_nan = pland_19_frac;
% pland_19_frac_nan(pland_19_frac_nan < 0.5) = NaN;
% pland_19_frac_nan(pland_19_frac_nan > 0) = 1;
% CM_19_crop = [];
% for i = 1:((size(CM_19,3)))
%     CM_19_nan = CM_19(:,:,i).*pland_19_frac_nan;
%     CM_19_crop(:,:,i) = CM_19_nan;
% end

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
dd19 = dd19.*tmp2; % multiply data by distance to get a matrix

%% Control Lateral Mobility
%look through 5mm radial distances and determine lateral mobility for each
%band

% control first
rad_dist = 0:50:2000; % distance we will calculate for (from 0 to 3100 mm, by 5 mm).. we could change this if we want finer grained information
% initialize empty matrices 
latmob90_18 = NaN(length(rad_dist)-1,1); %this will be mean hour it takes channel to visit 90% of radial band 
mean_latmob90_18 = NaN(length(rad_dist)-1,1); % mean 90% lateral mobility
std_latmob90_18 = NaN(length(rad_dist)-1,1); % standard 90% deviation mobility
latmob50_18 = NaN(length(rad_dist)-1,1); %this will be mean hour it takes channel to visit 50% of radial band
mean_latmob50_18 = NaN(length(rad_dist)-1,1); % mean 50% lateral mobility
std_latmob50_18 = NaN(length(rad_dist)-1,1); % standard deviation 50% lateral mobility

% radial loop
for k = 1:length(rad_dist)-1%loop to run through different radial distances from the end of the entrance channel.
    k
    %section to find x,y nodes for radial transect and generate matrix of
    %cross section topo and channel (yes/no) data
    radial_dd = (dd18 >= rad_dist(k) & dd18 < rad_dist(k+1));
    radial_dd = radial_dd.*basin18; % to remove radial transect outside basin
    area = sum(radial_dd(:), 'omitnan')*2.5*10^-5; %area in m2
    lat_chan_mob_18 = [];
    for i = 1:((size(CM_18,3)-60))
       for ii = i:size(CM_18,3)
           if ii == i
               CM_n = CM_18(:,:,ii).*radial_dd;
               lat_chan_mob_18_i = 1-((sum(CM_n(:), 'omitnan')*2.5*10^-5)/area);
               hour_18_i = (ii*2)-(i*2);
           else
               CM_n = CM_n + (CM_18(:,:,ii).*radial_dd);
               CM_n(CM_n > 1) = 1;
               lat_chan_mob_18_i = 1-((sum(CM_n(:), 'omitnan')*2.5*10^-5)/area);
               hour_18_i = ii-i;
               lat_chan_mob_18_i = [hour_18_i, lat_chan_mob_18_i];
               lat_chan_mob_18 = [lat_chan_mob_18;lat_chan_mob_18_i];
           end
        end
    end
    %for each timestep (1 hr, 2 hr, etc) calculate mean fraction not covered by flow
    [ud,ix,iy] = unique(lat_chan_mob_18(:,1));  
    mean_lat_chan_mob_18 = [ud, accumarray(iy,lat_chan_mob_18(:,2),[],@mean)];
    std_lat_chan_mob_18 = [ud, accumarray(iy,lat_chan_mob_18(:,2),[],@std)];
    
    %90% of area visited by channel
    [mean_hr] = find(mean_lat_chan_mob_18(:,2)<=0.1);
    if isempty(mean_hr)
        latmob90_18(k) = NaN;
        mean_latmob90_18(k) = NaN;
        std_latmob90_18(k) = NaN;
    else 
        hr = mean_hr(1);
        latmob90_18(k) = hr;
        mean_latmob90_18(k) = mean_lat_chan_mob_18(hr,2);
        std_latmob90_18(k) = std_lat_chan_mob_18(hr,2);
    end
    
    %50% of area visited by channel
    [mean_hr] = find(mean_lat_chan_mob_18(:,2)<=0.5);
    if isempty(mean_hr)
        latmob50_18(k) = NaN;
        mean_latmob50_18(k) = NaN;
        std_latmob50_18(k) = NaN;
    else 
        hr = mean_hr(1);
        latmob50_18(k) = hr(1);
        mean_latmob50_18(k) = mean_lat_chan_mob_18(hr,2);
        std_latmob50_18(k) = std_lat_chan_mob_18(hr,2);
    end
end

% Treatment Lateral Mobility
% initialize empty matrices 
rad_dist19 = 0:50:3000; % distance we will calculate for (from 0 to 3100 mm, by 5 mm).. we could change this if we want finer grained information
latmob90_19 = NaN(length(rad_dist19)-1,1); %this will be mean hour it takes channel to visit 90% of radial band 
mean_latmob90_19 = NaN(length(rad_dist19)-1,1); % mean 90% lateral mobility
std_latmob90_19 = NaN(length(rad_dist19)-1,1); % standard 90% deviation mobility
latmob50_19 = NaN(length(rad_dist19)-1,1); %this will be mean hour it takes channel to visit 50% of radial band
mean_latmob50_19 = NaN(length(rad_dist19)-1,1); % mean 50% lateral mobility
std_latmob50_19 = NaN(length(rad_dist19)-1,1); % standard deviation 50% lateral mobility

% radial loop
for k = 41:length(rad_dist19)-1%loop to run through different radial distances from the end of the entrance channel.
    k
    %section to find x,y nodes for radial transect and generate matrix of
    %cross section topo and channel (yes/no) data
    radial_dd = (dd19 >= rad_dist19(k) & dd19 < rad_dist19(k+1));
    radial_dd = radial_dd.*basin19;
    area = sum(radial_dd(:), 'omitnan')*2.5*10^-5; %area in m2
    lat_chan_mob_19 = [];
    for i = 1:((size(CM_19,3)-60))
       for ii = i:size(CM_19,3)
           if ii == i
               CM_n = CM_19(:,:,ii).*radial_dd;
               lat_chan_mob_19_i = 1-((sum(CM_n(:), 'omitnan')*2.5*10^-5)/area);
               hour_19_i = (ii*2)-(i*2);
           else
               CM_n = CM_n + (CM_19(:,:,ii).*radial_dd);
               CM_n(CM_n > 1) = 1;
               lat_chan_mob_19_i = 1-((sum(CM_n(:), 'omitnan')*2.5*10^-5)/area);
               hour_19_i = ii-i;
               lat_chan_mob_19_i = [hour_19_i, lat_chan_mob_19_i];
               lat_chan_mob_19 = [lat_chan_mob_19;lat_chan_mob_19_i];
           end
        end
    end
    %for each timestep (1 hr, 2 hr, etc) calculate mean fraction not covered by flow
    [ud,ix,iy] = unique(lat_chan_mob_19(:,1));  
    mean_lat_chan_mob_19 = [ud, accumarray(iy,lat_chan_mob_19(:,2),[],@mean)];
    std_lat_chan_mob_19 = [ud, accumarray(iy,lat_chan_mob_19(:,2),[],@std)];
    
    %90% of area visited by channel
    [mean_hr] = find(mean_lat_chan_mob_19(:,2)<=0.1);
    if isempty(mean_hr)
        latmob90_19(k) = NaN;
        mean_latmob90_19(k) = NaN;
        std_latmob90_19(k) = NaN;
    else 
        hr = mean_hr(1);
        latmob90_19(k) = hr;
        mean_latmob90_19(k) = mean_lat_chan_mob_19(hr,2);
        std_latmob90_19(k) = std_lat_chan_mob_19(hr,2);
    end
    
    %50% of area visited by channel
    [mean_hr] = find(mean_lat_chan_mob_19(:,2)<=0.5);
    if isempty(mean_hr)
        latmob50_19(k) = NaN;
        mean_latmob50_19(k) = NaN;
        std_latmob50_19(k) = NaN;
    else 
        hr = mean_hr(1);
        latmob50_19(k) = hr(1);
        mean_latmob50_19(k) = mean_lat_chan_mob_19(hr,2);
        std_latmob50_19(k) = std_lat_chan_mob_19(hr,2);
    end
end

distance = rad_dist/1000; % m
distance = distance(2:end);
distance19 = rad_dist19/1000;
distance19 = distance19(2:end);
%%Lets plot the data now
fig = figure;
plot(distance', latmob50_18, 'b-')
hold on
plot(distance19', latmob50_19, 'g-')
xlim([0 2])
ylim([0 75])
legend('control', 'treatment')
xlabel('distance from apex (m)')
ylabel('lateral mobility (hrs)')
set(gca, 'XMinorTick', 'On', 'YMinorTick', 'On')
saveas(fig, '../figures/esurf_Figure7c.pdf')

fig = figure;
plot(distance', latmob90_18, 'b-')
hold on
plot(distance19', latmob90_19, 'g-')
xlim([0 2])
ylim([0 200])
legend('control', 'treatment')
xlabel('distance from apex (m)')
ylabel('lateral mobility (hrs)')
set(gca, 'XMinorTick', 'On', 'YMinorTick', 'On')
saveas(fig, '../figures/esurf_Figure7d.pdf')

tmptreat = latmob90_19(1:40,:);
fig = figure;
plot(distance', tmptreat./latmob90_18, 'b-')
