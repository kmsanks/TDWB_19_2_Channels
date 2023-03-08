clear all; close all
%% The purpose of this script is to calculate hydrodynamic backwater length (Lb)
% These data are plotted in Figure 4 and Table 1 in Sanks et al.
% (2023) submitted to Earth Surface Dynamics

% We will use both LiDAR data and channel maps for the below calculation
%% First, lets set our relative path
cd '../' %We should have opened this code from 'C:\Users\XXXXXXX\Documents\GitHub\TDWB_19_2_Channels\code'

%% Load data
% control
cd './data'
load('ZD_18.mat'); %topography; elevation data (mm)
load('CM_18.mat'); %channel maps
% treatment
load('ZD_19.mat'); %topography; elevation data (mm)
load('CM_19.mat'); %channel maps
cd '../code'

%% Set parameters
% control
nx_18 = size(ZD_18, 1); %number of x locations on map
ny_18 = size(ZD_18,2); %number of y locations on map
nt_18 = size(ZD_18,3); %number of time steps in data set
dt_18 = 1; %delta t of time steps (hr)
xentrance_18 = 109; %x grid node location of the entrance channel
yentrance_18 = 271; %y grid node location of the entrance channel

% treatment
CM_19 = CM_19(:,:,2:2:end); %since we only have z for every other hour, crop channel maps
nx_19 = size(ZD_19,1); %number of x locations on map
ny_19 = size(ZD_19,2); %number of y locations on map
nt_19 = size(ZD_19,3); %number of time steps in data set
dt_19 = 2; %delta t of time steps (hr)
xentrance_19 = 214; %x grid node location of the entrance channel (x is down dip)
yentrance_19 = 397; %y grid node location of the entrance channel (y is strike)

% both experiments 
baselevel_rr = 0.25; %base level rise rate (mm/hr)
ocean_zero = 25; %ocean elevation at beginning of experiment (mm)
dx = 5; %5 mm grid cells

%% First, we need to crop topography to the channels
% Make elevation screen, so flow can be clipped to the areas above sea level
% control
z18 = []; %intialize empty matrix
for i= 1:nt_18; % loop through time steps
    %What is sea level at time i
    sl = (i*baselevel_rr*dt_18)+25; %first dry scan is at hour 1, so i = i
    elevationmask = ZD_18(:,:,i);  % elevation at time i
    elevationmask(elevationmask == 0) = NaN; %set everything outside basin = NaN
    elevationmask_rslr = elevationmask - sl; %subtract sea level from elevation to get elevation relative to sea level
    z18(:,:,i) = elevationmask_rslr; %save new matrix
end
% treatment
z19 = []; %intitialize empty matrix
for i= 1:nt_19; % loop through time steps
    %What is sea level at time i
    sl = ((i-1)*baselevel_rr*dt_19)+25; %first dry scan is hour 0, so we need i - 1
    elevationmask = ZD_19(:,:,i); % elevation at time i
    elevationmask(elevationmask == 0) = NaN; %set everything below sea level = NaN
    elevationmask_rslr = elevationmask - sl; %subtract sea level from elevation to get elevation relative to sea level
    z19(:,:,i) = elevationmask_rslr; %save new matrix
end


% Remove area outside channels
% control
zc18 = []; %initialize empty matrix
for i =1:nt_18 % loop through time steps
  c18 = (CM_18(:,:,i)).*(z18(:,:,i)); %multiply channel and elevations to get elevation inside channel
  c18(c18 == 0) = NaN; %0 is outside channels 
  zc18(:,:,i) = c18; %save new matrix
end
% treatment
zc19 = []; %initialize empty matrix
for i =1:(nt_19-1) % loop through timesteps
  c19 = (CM_19(:,:,i)).*(z19(:,:,i)); %multiply channel and elevations to get elevation inside channel
  c19(c19 == 0) = NaN; %0 is outside channels 
  zc19(:,:,i) = c19; %save new matrix
end

%% Below, we will calculate elevation of the channel bed through time and associated backwater length

% First, we need distance to each pixel
% control
[X18, Y18] = meshgrid(1:ny_18, 1:nx_18); % get pixels
dd18 = sqrt((X18 - yentrance_18).^2 + (Y18 - xentrance_18).^2)*dx; % calculate distance to each pixel from apex

%make everything outside of basin a NaN
tmp = zeros(796,522); % empty matrix
z = ZD_18(:,:,1); %fill with elevations
z(z == 0.) = NaN; % everything outside basin is NaN
tmp2 = z.*tmp; %inside basin = 0, outside basin = NaN
tmp2(tmp2 == 0.) = 1; %inside basin = 1
dd18 = dd18.*tmp2; % get distance only inside the basin

% treatment
[X19, Y19] = meshgrid(1:ny_19, 1:nx_19); % get pixels
dd19 = sqrt((X19 - yentrance_19).^2 + (Y19 - xentrance_19).^2)*dx; % calculate distance to each pixel from apex

%make everything outside of basin a NaN
tmp = zeros(750,747); % empty matrix
z = ZD_19(:,:,1); %fill with elevations
tmp2 = z.*tmp; % inside basin = 0, outside basin = NaN
tmp2(tmp2 == 0.) = 1; %inside basin = 1
dd19 = dd19.*tmp2; % get distance only inside the basin

% Now we are going to loop through to find distance to backwater
backwater_elev18 = []; %intialize empty matrix
backwater_point18 = []; %initialize empty matrix
backwater_end18 = []; %initialize empty matrix
dist = 0:10:3100; % distances to loop through from 0 to 3100 mm every 10 mm
backwater_mat_18 = NaN(length(dist)-1, 2, size(CM_18,3)); %matrix to input values
for i = 1:size(CM_18,3) % loop through every hour
    backwater = NaN(length(dist)-1, 3);
    for k = 1:(length(dist)-1) %loop through distances
        idx = dd18(dd18 >= dist(k) & dd18 < dist(k+1)); 
        radial_dd = dd18 >= dist(k) & dd18 < dist(k+1); %10 mm radial transect of interest
        is_shot = CM_18(:,:,i).*radial_dd; %pixels in transect that are in the channel
        is_shot(is_shot == 0.) = NaN; %everything else NaN
        z = zc18(:,:,i).*is_shot; %elevation relative to sea level in channel in transect
        %now we want to determine if at least 16% of the elevations are at
        %or below sea level
        el = prctile(z(:),16); % 16 percentile elevation of channel transect
   % This calculated first elevation below sea level and distance
        %if el <= 0 %transect is below sea level
        %   backwater_elev = [backwater_elev; el]; %save elevation of transect
        %   backwater_dist = [backwater_dist; dist(k)]; %save radial distance of transect
        %end
   % We actually want to find the max consecutive distance
        backwater(k,1) = el; %save elevation of transect to tmp variable
        backwater(k,2) = dist(k); %save radial distance of transect to tmp variable
        if el>0
            backwater(k,3) = 0;
        else
            backwater(k,3) = 1;
        end
        backwater_mat_18(k,1,i) = dist(k); %save to matrix
        backwater_mat_18(k,2,i) = el; %save to matrix
    end
  % This saved that first elevation below sea level and distance  
    %if isempty(backwater_elev)
    %   backwater_elev18(i) = NaN; %save elevation = NaN if channel bed is not below sea level
    %   backwater_dist18(i) = NaN; %save distance = NaN if channel bed is not below sea level
    %else
    %   backwater_elev18(i) = backwater_elev(1); %save elevation if bed is below sea level
    %   backwater_dist18(i) = backwater_dist(1); %save distance if bed is below sea level
    %end 
 % We actually want to find max consecutive distance from end of channel
   backwater(any(isnan(backwater), 2), :) = [];
   n = backwater(:,3); %logical of channel bed elevation above/below sea level
   props = regionprops(logical(n),'Area','PixelIdxList'); %find the number of 1s and the indicies
   bw_sect = NaN(length(props), 3); %initialize matrix
   for  k = 1:length(props) %loop through groups of consectuive 1s
        bw_sect(k,1) = props(k).Area(1); %save number of consecutive 1s for each group
        bw_sect(k,2) = props(k).PixelIdxList(1); %save start index
        bw_sect(k,3) = props(k).PixelIdxList(end);%save end index
   end
   %Here we will find maximum cosecutive backwater length
   if isempty(bw_sect) %if channel bed is never below sea level
   %if  backwater(size(backwater,1),1)>=0
       backwater_elev18(i) = NaN; %save elevation = NaN if channel tip is not below sea level
       backwater_point18(i) = NaN; %save distance = NaN if channel tip is not below sea level
       backwater_end18(i) = NaN; %save distance = NaN if channel tip is not below sea level
   else %if channel bed does go below sea level
        val_max = max(bw_sect, [], 1); %[conseutive ones, index start, index end] for max consecutive ones
        backwater_elev18(i) = backwater(val_max(2),1); %beginning elevation
        backwater_point18(i) = backwater(val_max(2),2); %backwater point
        backwater_end18(i) = backwater(val_max(3),2); %end distance
   end
end

% treatment
backwater_elev19 = []; %initialize empty matrix
backwater_point19 = []; %initialize empty matrix
backwater_end19 = []; %initialize empty matrix
backwater_mat_19 = NaN(length(dist)-1, 2, size(CM_19,3)); %matrix to input values
for i = 1:size(CM_19,3) %loop through every hour
    backwater = NaN(length(dist)-1, 3);
    for k = 1:(length(dist)-1) %loop through distances
        idx = dd19(dd19 >= dist(k) & dd19 < dist(k+1)); %10 mm radial transect of interest
        radial_dd = dd19 >= dist(k) & dd19 < dist(k+1); %pixels in transect that are in the channel
        is_shot = CM_19(:,:,i).*radial_dd; %pixels in transect that are in the channel
        [row, col] = find(~isnan(is_shot)); %find time steps where there is no channel map
        if isempty(row) %if no channel map, fill array with NaN 
            backwater_elev19(i) = NaN; 
            backwater_point19(i) = NaN;
            backwater_end19(i) = NaN;
            backwater_length19(i) = NaN;
            backwater_mat_19(k,1,i) = NaN;
            backwater_mat_19(k,2,i) = NaN;
        else %otherwise find elevations
            is_shot(is_shot == 0.) = NaN;  %everything outside of channel in transect is NaN
            z = z19(:,:,i).*is_shot; %elevation relative to sea level
            %now we want to determine if at least 16% of the elevations are at
            %or below sea level
            el = prctile(z(:),16);
            %if el <= 0  %transect is below sea level
            %   backwater = [backwater; el]; %save elevation of transect
            %   backwater_dist = [backwater_dist; dist(k)]; %save radial distance of transect
            %end
            backwater(k,1) = el; %save elevation of transect to tmp variable
            backwater(k,2) = dist(k); %save radial distance of transect to tmp variable
            if  el>0
                backwater(k,3) = 0;
            else
                backwater(k,3) = 1;
            end
            backwater_mat_19(k,1,i) = dist(k); %save to matrix
            backwater_mat_19(k,2,i) = el; %save to matrix
        end
    end
    % We actually want to find max consecutive distance from end of channel
    backwater(any(isnan(backwater), 2), :) = [];
    n = backwater(:,3); %logical of channel bed elevation above/below sea level
    props = regionprops(logical(n),'Area','PixelIdxList'); %find the number of 1s and the indicies
    bw_sect = NaN(length(props), 3); %initialize matrix
    for  k = 1:length(props) %loop through groups of consectuive 1s
         bw_sect(k,1) = props(k).Area(1); %save number of consecutive 1s for each group
         bw_sect(k,2) = props(k).PixelIdxList(1); %save start index
         bw_sect(k,3) = props(k).PixelIdxList(end);%save end index
    end
    %Here, we will find maximum consecutive backwater length
    if   isempty(bw_sect) %if channel bed is never below sea level
    %if  backwater(size(backwater,1),1)>=0
         backwater_elev19(i) = NaN; %save elevation = NaN if channel tip is not below sea level
         backwater_point19(i) = NaN; %save distance = NaN if channel tip is not below sea level
         backwater_end19(i) = NaN; %save distance = NaN if channel tip is not below sea level
    else %if channel bed does go below sea level
         val_max = max(bw_sect, [], 1); %[conseutive ones, index start, index end] for max consecutive ones
         backwater_elev19(i) = backwater(val_max(2),1); %beginning elevation
         backwater_point19(i) = backwater(val_max(2),2); %backwater point
         backwater_end19(i) = backwater(val_max(3),2); %end distance
    end
end

%% Backwater length
% control
% calculate backwater length
% 0 if backwater is only one pixel, NaN if no backwater
% because of this we will add 10 mm (width of the transect) to all of the calculated values
backwater_length18 = (backwater_end18-backwater_point18)+10;
backwater_length18(isnan(backwater_length18))=0; %replace no backwater with NaN

% treatment
% 0 if backwater is only one pixel, NaN if no backwater or no channel map
% because of this we will add 10 mm (width of the transect) to all of the calculated values
backwater_length19 = (backwater_end19-backwater_point19)+10; 
% Find timesteps with no channel maps in treatment to plot as "x's" on plot and remove from statistical calculations
% which timesteps have no channel maps so we can plot later
% replace backwater length with NaN
hour_nomaps19 = [];
bw_nomaps19 = [];
for i = 1:size(CM_19,3)
    maps = CM_19(:,:,i);
    if sum(maps(:), 'omitnan') == 0
        hour = i;
        hour_nomaps19 = [hour_nomaps19;hour];
        bw_nomaps19 = [bw_nomaps19;0];
    end
end

% check number of backwater NaN and number of no maps
isequal(length(bw_nomaps19),sum(isnan(backwater_length19)))
% they are equal, so all timesteps have backwater besides ones without
% channel maps: no need to replace any data in the array
%% Calculate some statistics:
%backwater stats to report in manuscript
frac_backwater_control = sum(backwater_length18>0)/nt_18; %fraction of time backwater
greater_500_control = sum(backwater_length18>500)/nt_18; %fraction of time backwater is > 500 mm
Lb_mean_control = mean(backwater_length18)/1000; %m
Lb_std_control = std(backwater_length18)/1000; %m

t19 = 280-length(bw_nomaps19);
frac_nomaps_treat = sum(isnan(backwater_length19))/280; %fraction of time we are missing channel maps
frac_backwater_treat = sum(backwater_length19>0)/t19; %fraction of time there is backwater
greater_500_treat = sum(backwater_length19>500)/t19; %fraction of time backwater is > 0.5
Lb_mean_treat = mean(backwater_length19, 'omitnan')/1000; %m
Lb_std_treat = std(backwater_length19, 'omitnan')/1000; %m

%% Now, we want to plot the data
% plot through nans for treatment data
y = backwater_length19/1000;
x = 1:2:nt_18;
idx_plot_19 = ~any(isnan(y),1);

fig1 = figure();
plot(1:nt_18, (backwater_length18/1000), 'b-', 'LineWidth', 2)
hold on
plot(x(idx_plot_19),y(idx_plot_19), 'g-', 'LineWidth', 2)
yline((mean(backwater_length18, 'omitnan')/1000), 'b:', 'LineWidth', 2)
yline((mean(backwater_length19, 'omitnan')/1000), 'g:', 'LineWidth', 2)
plot(hour_nomaps19*2, bw_nomaps19, 'kx', 'LineWidth', 1.5, 'MarkerSize', 10)
xlim([0 560])
xlabel('time (hrs)')
ylabel('backwater length (m)')
legend('control', 'treatment', 'control mean', 'treatment mean', 'no channel map treatment')
set(gca,'XMinorTick','on','YMinorTick','on')
y_width=7.25;x_width=9.125
set(gcf, 'PaperPosition', [0 0 x_width y_width]);
saveas(fig1, '../figures/esurf_Figure4.pdf')
