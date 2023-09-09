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
load('ZD_18.mat') % control topography (mm)
load('ZD_19.mat') % treatment topography from Sam (mm)
load('CM_18.mat') % control channel maps (binary, 1s = channels)
load('CM_19.mat') % treatment channel maps (binary, 1s = channels)
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
CM_18 = CM_18(:,:,1:2:560); % hour 1, 3, etc.
ZD_19 = ZD_19(:,:,2:end); % remove first time step because that is = 0 and we want to start at t = 1
CM_19 = CM_19(:,:,2:2:end); % remove first time step 

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

%% We only want to calculate for area above -9mm relative to sea level, so let's mask the basin
% control
% terr_area18 = []; % initialize terrestrial area binary matrix
% for i= 1:nt_18 % loop through all timesteps
%     %What is sea level at time i
%     sl = ((i-1)*baselevel_rr)+ocean_zero; % first scan is hour 1, so sl = 25.25 mm
%     elevationmask = ZD_18(:,:,i); 
%     elevationmask(elevationmask == 0) = NaN;
%     elevationmask_rslr = elevationmask - sl;
%     elevationmask_rslr(elevationmask_rslr < -9) = NaN;
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
%     terr_area18(:,:,i) = createMask(terr);
% end
% 
% % treatment
% terr_area19 = []; % initialize terrestrial area binary matrix
% for i= 1:nt_19 % loop through all timesteps
%     %What is sea level at time i
%     sl = ((i-1)*baselevel_rr)+ocean_zero; % first scan is hour 1, so sl = 25.25 mm
%     elevationmask = ZD_19(:,:,i); 
%     elevationmask(elevationmask == 0) = NaN;
%     elevationmask_rslr = elevationmask - sl;
%     elevationmask_rslr(elevationmask_rslr < -9) = NaN;
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
%     terr_area19(:,:,i) = createMask(terr);
% end
% 
% save('../data/terr_area18.mat', 'terr_area18') % terrestrial area binary > -9mm rsl
% save('../data/terr_area19.mat', 'terr_area19') % terrestrial area binary > -9mm rsl

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
zrsl18 = [];
for i= 1:nt_18
    %What is sea level at time i
    sl = ((i-1)*baselevel_rr)+ocean_zero; %first wet scan starts at hour 1
    elevationmask = ZD_18(:,:,i); 
    elevationmask(elevationmask == 0) = NaN;
    elevationmask_rslr = elevationmask - sl;
    elevationmask_rslr(elevationmask_rslr < -50) = NaN;
    zrsl18(:,:,i) = elevationmask_rslr;
end

% treatment
zrsl19 = [];
for i= 1:nt_19
    %What is sea level at time i
    sl = ((i-1)*baselevel_rr)+ocean_zero; %first scan is t = 2, which is hour 1
    elevationmask = ZD_19(:,:,i); 
    elevationmask(elevationmask == 0) = NaN;
    elevationmask_rslr = elevationmask - sl;
    elevationmask_rslr(elevationmask_rslr < -50) = NaN;
    zrsl19(:,:,i) = elevationmask_rslr;
end

%% Channel properties loop: here we will calculate radial channel properties (channel area, width, number of channels)
% control first
rad_dist = 0:50:3100; % distance we will calculate for (from 0 to 3100 mm, by 50 mm or 5 cm).. we could change this if we want finer grained information
% initialize empty matrices 
width18 = NaN(length(rad_dist),nt_18); % trunk channel width
n_chan18 = NaN(length(rad_dist),nt_18); % number of channels
meanelev18 = NaN(length(rad_dist),1); % mean channel bed elevation relative to sea level
stdelev18 = NaN(length(rad_dist),1); % standard deviation of channel bed elevation relative to sea level

% loop through radial distances 
for k = 1:(length(rad_dist)-1) % loop to run through different radial distances from the apex.
    % print which radial transect the loop is on
    caption = sprintf('Radial segment %d', k);
    caption2 = sprintf(' of %d', length(rad_dist)-1);
    fprintf('%s\n', strcat(caption,caption2))
    % section to find x,y nodes for radial transect and generate matrix of
    % cross section topo and channel (yes/no) data
    radial_dd = dd18 >= rad_dist(k) & dd18 < rad_dist(k+1); % can look at this using imagesc(radial_dd) to visualize a 0.1 m radial transect
    
    bed = [];
    % loop through time
    for i = 1:size(CM_18,3) % loop through all timesteps
        area = terr_area18(:,:,i).*basin18;
        area(isnan(area))=0; 
        is_shot = CM_18(:,:,i).*radial_dd.*area; % in channel or no? on delta >-9 mm or no?
        is_chan = is_shot; % need 0s and 1s for bwlabel
        is_shot(is_shot == 0) = NaN; % need NaNs for z_rsl
        z_rsl = zrsl18(:,:,i).*is_shot; % elevation of channel relative to sea level in radial transect
        [label,n] = bwlabel(is_chan); % label gives a unique number to each individual segement recognized; n is number of channels (or segments)
        rad = area.*radial_dd; % radial transect
        radial_area = sum(rad(:),'omitnan')*0.25; % area of radial transect in cm^2
        % save total channel area for each distance through time
        ca = [];
        width = [];
        % loop through channel segments
        for j = 1:n
            % binary channel segments
            tmp = label;
            tmp(tmp ~= j) = NaN; % remove all data that is not in segement n
            tmp(tmp > 0) = 1; % turn segement data to 1
            tmp(tmp < 1) = NaN; % make everything else NaN
            
            % channel area
            carea = sum(tmp(:), 'omitnan')*0.25; % area of channel segment in cm^2
            ca = [ca,carea]; % save area for each segment
            
            % channel width (take out channel tips?)
            w = carea/5; %width in cm assuming length is 5 cm (length of radial bin)
            width = [width,w]; % save channel width for each segment

            % channel depth
            z_tmp = z_rsl.*tmp; % segment elevations to get depth
            tmp_bed = prctile(z_tmp(:), 5); % 5th% channel elevation for each segment
            % save channel bed elevation for each distance and each channel through time
            bed = [bed,tmp_bed]; % to save all channel depths for each radial transect
        end 
   
        % save channel area and fraction
        if isempty(ca)
           width18(k,i) = NaN;
        else
           width18(k,i) = max(width);
        end    
        
        % save number of channels for each distance through time
        if isempty(n)
           n_chan18(k,i) = NaN;
        else
           n_chan18(k,i) = n;
        end
    end 
    % save channel bed elevations
    meanelev18(k) = mean(bed, 'omitnan');
    stdelev18(k) = std(bed, 'omitnan');
end 

% treatment
% replace channel maps for treatment with NaN if time step does not have one
for i = 1:size(CM_19,3) 
    if sum(sum(CM_19(:,:,i), 'omitnan'), 'omitnan') == 0
        CM_19(:,:,i) = NaN;
    end
end

% initialize empty matrices 
width19 = NaN(length(rad_dist),nt_19); % trunk channel width
n_chan19 = NaN(length(rad_dist),nt_19); % number of channels
meanelev19 = NaN(length(rad_dist),1); % mean channel bed elevation relative to sea level
stdelev19 = NaN(length(rad_dist),1); % standard deviation of channel bed elevation relative to sea level
% run loop
for k = 1:(length(rad_dist)-1)%loop to run through different radial distances from the end of the apex.
    % print which radial segment the loop is on
    caption = sprintf('Radial segment %d', k);
    caption2 = sprintf(' of %d', length(rad_dist)-1);
    fprintf('%s\n', strcat(caption,caption2))

    %section to find x,y nodes for radial transect and generate matrix of
    %cross section topo and channel (yes/no) data
    idx = dd19(dd19 >= rad_dist(k) & dd19 < rad_dist(k+1));
    radial_dd = dd19 >= rad_dist(k) & dd19 < rad_dist(k+1);
    
    bed = [];
    % loop through time
    for i = 1:size(CM_19,3) % loop through all timesteps
        area = terr_area19(:,:,i).*basin19;
        area(isnan(area))=0; 
        is_shot = CM_19(:,:,i).*radial_dd.*area; % in channel or no? on delta >-9 mm or no?
        is_chan = is_shot; % need 0s and 1s for bwlabel

        % is there a channel map?
        if sum(is_shot(:), 'omitnan') == 0 % do not calculate channel properties for timesteps with no channel map
            area19(k,i) = NaN;
            chanfrac19(k,i) = NaN;
            width19(k,i) = NaN;

        else % calculate channel properties
            is_shot(is_shot == 0) = NaN; % need NaNs for z_rsl
            z_rsl = zrsl19(:,:,i).*is_shot; % elevation of channel relative to sea level in radial transect
            [label,n] = bwlabel(is_chan); % label gives a unique number to each individual segement recognized; n is number of channels (or segments)
            rad = area.*radial_dd; % radial transect
            radial_area = sum(rad(:),'omitnan')*0.25; % area of radial transect in cm^2
            
            % loop through channel segments
            ca = [];
            width = [];
            for j = 1:n
                % binary channel segments
                tmp = label;
                tmp(tmp ~= j) = NaN; % remove all data that is not in segement n
                tmp(tmp > 0) = 1; % turn segement data to 1
                tmp(tmp < 1) = NaN; % make everything else NaN
                
                % channel area
                carea = sum(tmp(:), 'omitnan')*0.25; % area of channel segment in cm^2
                ca = [ca,carea]; % save area for each segment
                
                % channel width
                w = carea/5; %width in cm assuming length is 5 cm (length of radial bin)
                width = [width,w]; % save channel width for each segment
    
                % channel depth
                z_tmp = z_rsl.*tmp; % segment elevations to get depth
                tmp_bed = prctile(z_tmp(:), 5); % 5th% channel elevation for each segment
                % save channel bed elevation for each distance and each channel through time
                bed = [bed,tmp_bed]; % to save all channel depths for each radial transect
            end 
   
            % save channel area and fraction
            if isempty(ca)
               width19(k,i) = NaN;
            else
               width19(k,i) = max(width);
            end    
            
            % save number of channels for each distance through time
            if isempty(n)
               n_chan19(k,i) = NaN;
            else
               n_chan19(k,i) = n;
            end
        end 
    end 
    % save channel bed elevations
    meanelev19(k) = mean(bed, 'omitnan');
    stdelev19(k) = std(bed, 'omitnan');   
end 

%% Calculate statistics for Table 1:
% trunk channel width
width_mean18 = mean(width18(:), 'omitnan');
width_std18 = std(width18(:), 'omitnan');
width_mean19 = mean(width19(:), 'omitnan');
width_std19 = std(width19(:), 'omitnan');

%% Plot the data: Figure 4a
ybars = [-9 5]; % marsh window
rad_dist = 0:0.05:3.1;

array18 = [rad_dist; meanelev18'; stdelev18'];
cols = any(isnan(array18),1);
array18(:,cols) = [];

array19 = [rad_dist; meanelev19'; stdelev19'];
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
alpha(0.1)
yline(0, 'k-', 'linewidth', 2)
plot(x18, y18, 'b', 'LineWidth', 2)
plot(x19, y19, 'g', 'LineWidth', 2)
ylim([-15 35])
xlim([0 3])
ylabel('mean channel bed elevation relative to sea level (mm)')
xlabel('distance from apex (m)') 
legend('control mean', 'treatment mean', 'control stdev', 'treatment stdev', 'marsh window', 'sea level')
set(gca, 'XMinorTick', 'On', 'YMinorTick', 'On')
set(gcf, 'PaperUnits', 'inches');
y_width=7.25;x_width=9.125;
set(gcf, 'PaperPosition', [0 0 x_width y_width]);
saveas(fig, '../figures/esurf_Figure4a.pdf')

%% Plot the data: Figure 3c
% create arrays
w18 = mean(width18,2, 'omitnan');
stdev18 = std(width18,[],2, 'omitnan');
chan_array18 = [rad_dist; w18'; stdev18'];
cols = any(isnan(chan_array18),1);
chan_array18(:,cols) = [];

w19 = mean(width19,2, 'omitnan');
stdev19 = std(width19,[],2, 'omitnan');
chan_array19 = [rad_dist; w19'; stdev19'];
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
plot(x18, y18, 'b-', 'LineWidth', 2)
hold on
plot(x19, y19, 'g-', 'LineWidth', 2)
patch([x18 fliplr(x18)], [curve1_18 fliplr(curve2_18)], 'b')
patch([x19 fliplr(x19)], [curve1_19 fliplr(curve2_19)], 'g')
alpha(0.10)
plot(x18, y18, 'b-', 'LineWidth', 2)
plot(x19, y19, 'g-', 'LineWidth', 2)
ylabel('max channel width (cm)')
ylim([0 18])
xlim([0 3])
xlabel('distance from apex (m)') 
legend('control mean', 'treatment mean', 'control stdev', 'treatment stdev')
set(gca, 'YMinorTick', 'On', 'XMinorTick', 'On')
set(gcf, 'PaperUnits', 'inches');
y_width=7.25;x_width=9.125;
set(gcf, 'PaperPosition', [0 0 x_width y_width]);
saveas(fig, '../figures/esurf_Figure4c.pdf')

%% Plot the data: Figure 3d
% create arrays
n18 = mean(n_chan18,2, 'omitnan');
stdev18 = std(n_chan18,[],2, 'omitnan');
chan_array18 = [rad_dist; n18'; stdev18'];
cols = any(isnan(chan_array18),1);
chan_array18(:,cols) = [];

n19 = mean(n_chan19,2, 'omitnan');
stdev19 = std(n_chan19,[],2, 'omitnan');
chan_array19 = [rad_dist; n19'; stdev19'];
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
plot(x18, y18, 'b-', 'LineWidth', 2)
hold on
plot(x19, y19, 'g-', 'LineWidth', 2)
patch([x18 fliplr(x18)], [curve1_18 fliplr(curve2_18)], 'b')
patch([x19 fliplr(x19)], [curve1_19 fliplr(curve2_19)], 'g')
alpha(0.10)
plot(x18, y18, 'b-', 'LineWidth', 2)
plot(x19, y19, 'g-', 'LineWidth', 2)
ylim([0 5])
xlim([0 3])
ylabel('number of channels')
xlabel('distance from apex (m)') 
legend('control mean', 'treatment mean', 'control stdev', 'treatment stdev')
set(gca, 'YMinorTick', 'On', 'XMinorTick', 'On')
set(gcf, 'PaperUnits', 'inches');
y_width=7.25;x_width=9.125;
set(gcf, 'PaperPosition', [0 0 x_width y_width]);
saveas(fig, '../figures/esurf_Figure4d.pdf')

%% Now we will calculate channel depth, aggradation, and channel in-filling and compensation timescales
% Lets buffer the channels so levee sedimentation is included in channel aggradation 
se = strel('square',4); %create a square of 4 pixels (2 cm) around each 1 (channel) pixel, this will buffer our channel maps to ensure we get levee crests
% control
CM_buffer_18 = [];
for j = 1:nt_18 %go by 2 so no saddler effects here
    buffer = imdilate(CM_18(:,:,j), se); %buffer channel map by 10 pixels
    CM_buffer_18(:,:,j) = buffer;
end 
FF_buffer_18 = ~CM_buffer_18; % inverse of channel buffer
FF_buffer_18 = double(FF_buffer_18); % convert to double
area18 = (terr_area18(:,:,1:279) + CM_buffer_18(:,:,1:279)).*basin18;
area18(area18>1)=1; % channel or >-9mm rsl = 1
area18(area18==0)=NaN; % turn 0s back to NaN
CM_buffer_18(CM_buffer_18==0)=NaN; % turn outside channel to NaN
FF_buffer_18(FF_buffer_18==0)=NaN; % turn channels to NaN

% treatment
CM_buffer_19 = [];
for j = 1:nt_19
    buffer = imdilate(CM_19(:,:,j), se);
    CM_buffer_19(:,:,j) = buffer;
end 
FF_buffer_19 = ~CM_buffer_19; % inverse of channel buffer
FF_buffer_19 = double(FF_buffer_19); % convert to double
area19 = (terr_area19(:,:,1:279) + CM_buffer_19(:,:,1:279)).*basin19;
area19(area19>1)=1; % channel or >-9mm rsl = 1
area19(area19==0)=NaN; % turn 0s back to NaN
CM_buffer_19(CM_buffer_19==0)=NaN; % turn outside channel to NaN
FF_buffer_19(FF_buffer_19==0)=NaN; % turn channels to NaN

%% Clean treatment channel maps

% depth
% replace channel maps for treatment with NaN if time step does not have one
for i = 1:size(CM_buffer_19,3) 
    if sum(sum(CM_buffer_19(:,:,i), 'omitnan'), 'omitnan') == -Inf
        CM_buffer_19(:,:,i) = NaN;
        FF_buffer_19(:,:,i) = NaN;
        area19(:,:,i) = NaN;
    end
end

% % aggradation
% % Replace timesteps with no channel maps with the channel map from the
% % next time step for the treatment experiment
% for i = (nt_19-1):-1:1 %I know 280 has a channel map, so I can start at 279 and replace with channel map that comes after, this will work for ones that have multiple no maps in a row
%     if sum(sum(CM_buffer_19(:,:,i), 'omitnan'), 'omitnan') == 0
%         CM_buffer_19(:,:,i) = CM_buffer_19(:,:,i+1);
%     end
% end
%% Calculate channel and far-field aggradation
% Total dz: control
dZD_18 = diff(ZD_18,1,3); % difference the maps along time dimension (3)
dZD_18(dZD_18 > 45) = NaN; % remove extraneous aggradation values
dZD_18(dZD_18 < -15) = NaN; % remove extraneous erosion values
% Total dz: treatment
dZD_19 = diff(ZD_19,1,3); % difference the maps along time dimension (3)
dZD_19(dZD_19 > 45) = NaN; % remove extraneous aggradation values
dZD_19(dZD_19 < -15) = NaN; % remove extraneous erosion values

% When you difference the maps there is one less timestep
CM_buffer_18 = CM_buffer_18(:,:,1:279); 
FF_buffer_18 = FF_buffer_18(:,:,1:279);  
CM_buffer_19 = CM_buffer_19(:,:,1:279);  
FF_buffer_19 = FF_buffer_19(:,:,1:279);  

% Far field maps (remove channel from topo): control
dZD_18_FF_buffer = FF_buffer_18.*dZD_18.*area18; % difference elevation of far field

% Far field maps (remove channel from topo): treatment
dZD_19_FF_buffer = FF_buffer_19.*dZD_19.*area19; % difference elevation of far field

% Channel maps (remove far field from topo): control
dZD_18_chan_buffer = CM_buffer_18.*dZD_18.*area18; % elevation of channels

% Channel maps (remove far field from topo): treatment
dZD_19_chan_buffer = CM_buffer_19.*dZD_19.*area19; % elevation of channels

%% Calculate aggradation statistics for Table 1
mean_chan_agg18 = median(dZD_18_chan_buffer(:), 'omitnan')/2; %mm/hr
stdev_chan_agg18 = std(dZD_18_chan_buffer(:), 'omitnan')/2; %mm/hr

mean_chan_agg19 = median(dZD_19_chan_buffer(:), 'omitnan')/2; %mm/hr
stdev_chan_agg19 = std(dZD_19_chan_buffer(:), 'omitnan')/2; %mm/hr

mean_ff_agg18 = median(dZD_18_FF_buffer(:), 'omitnan')/2; %mm/hr
stdev_ff_agg18 = std(dZD_18_FF_buffer(:), 'omitnan')/2; %mm/hr

mean_ff_agg19 = median(dZD_19_FF_buffer(:), 'omitnan')/2; %mm/hr
stdev_ff_agg19 = std(dZD_19_FF_buffer(:), 'omitnan')/2; %mm/hr

%% Calculate aggradation, compensation timescale, and channel in-filling rate

% initialize empty matrices to fill 
rad_dist = 0:5:3100; % each radial transect is 50 mm or 5 cm
Hc_18 = NaN(length(rad_dist),nt_18-1); % channel depth with distance and time to use in aggradation below
Hc_tot18 = NaN(length(rad_dist),100); % to save percentile channel depth with distance
Hc_all18 = []; % all channel depths
zff_18 = NaN(length(rad_dist),nt_19-1); % to save mean far-field elevation rsl
zc_18 = NaN(length(rad_dist),nt_19-1); % to save channel bed elevation rsl
zlev_18 = NaN(length(rad_dist),nt_19-1); % to save levee crest elevation rsl
agg18_mean = NaN(length(rad_dist),nt_18-1); % total aggradation rate for distance through time
c_agg18_mean = NaN(length(rad_dist),nt_18-1); % channel aggradation rate for distance through time
ff_agg18_mean = NaN(length(rad_dist),nt_18-1); % far-field aggradation rate for distance through time
d_agg18_mean = NaN(length(rad_dist),nt_18-1); % difference in mean aggradation between channel and ff for distance through time
agg18_median = NaN(length(rad_dist),nt_18-1); % total aggradation rate for distance through time
c_agg18_median = NaN(length(rad_dist),nt_18-1); % channel aggradation rate for distance through time
ff_agg18_median = NaN(length(rad_dist),nt_18-1); % far-field aggradation rate for distance through time
d_agg18_median = NaN(length(rad_dist),nt_18-1); % difference in mean aggradation between channel and ff for distance through time
% loop through radial distances 
for k = 2:(length(rad_dist)-1) % loop to run through different radial distances from the end of the entrance channel.
    % print which radial segment the loop is on
    caption = sprintf('Radial segment %d', k);
    caption2 = sprintf(' of %d', length(rad_dist)-1);
    fprintf('%s\n', strcat(caption,caption2))

    % section to find x,y nodes for radial transect and generate matrix of
    % cross section topo and channel (yes/no) data
    idx = dd18(dd18 >= rad_dist(k) & dd18 < rad_dist(k+1));
    radial_dd = dd18 >= rad_dist(k) & dd18 < rad_dist(k+1); % can look at this using imagesc(radial_dd) to visualize a 0.1 m radial transect
    % initialize empty matrices 
    depth = [];
    for i = 1:(size(CM_buffer_18,3))
        % only want area on the delta top >-9 mm rsl or in the channels 
        CMtmp = CM_buffer_18(:,:,i); % we need tmp channel map
        CMtmp(isnan(CMtmp))=0; % to turn NaNs back to 0
        area = (terr_area18(:,:,i) + CMtmp).*basin18; %so we can add to area
        area(area>1)=1; % channel or >-9mm rsl = 1
        area(area==0)=NaN; % turn 0s back to NaN
        
        % dz
        rads = double(radial_dd); % turn into double
        rads(rads==0)=NaN; % so we can turn 0s to ones
        dz = dZD_18(:,:,i).*rads.*area; % so we can get dz only in the radius
        
        % channel and ff data in radial transect
        dz_chan = dZD_18_chan_buffer(:,:,i).*rads.*area; % multiply by radial transect 
        dz_ff = dZD_18_FF_buffer(:,:,i).*rads.*area; % multiply by radial transect
        
        % elevation data for channel depths
        z = zrsl18(:,:,i); % elevation data for channel depth
        is_shot = CM_buffer_18(:,:,i).*rads.*area; % in channel or no?
        z = z.*is_shot; % elevation in the channels
        is_shot(isnan(is_shot)) = 0; % can't have NaNs in bwlabel
        
        % mean far-field elevation relative to sea level
        ztmp = zrsl18(:,:,i).*rads.*area.*FF_buffer_18(:,:,i); % only in terrestrial area > -9mm or in channel
        if sum(~isnan(ztmp(:)))==0 % if there are no far-field pixels in transect
            zff_18(k,i) = NaN; % save NaN
        else
            zff_18(k,i) = median(ztmp(:), 'omitnan'); % otherwise save median far-field elevation (mm rsl)
        end

        % how many individual channel segments ?
        [label,n] = bwlabel(is_shot); % label gives a unique number to each individual segement recognized; n is number of channels (or segments)
       
        % channel depth
        % we will calculate depth for each channel here
        depthmax = [];
        zc = [];
        zlev = [];
        for j = 1:n % loop through channel segments
            tmp = label;
            tmp(tmp ~= j) = NaN; % remove all data that is not in segement n
            tmp(tmp > 0) = 1; % turn segement data to 1
            tmp(tmp < 1) = NaN; % make everything else NaN
            z_tmp = z.*tmp; % segment elevations to get depth
            zs = z_tmp(~isnan(z_tmp)); % remove data not in the channel
            tmp_depth = (max(zs) - min(zs)); % channel depth in mm; we will take the max elevation as levee elevation and min elevation as thalweg                  
            zc_tmp = min(zs); % channel bed elevation relative to sea level
            zlev_tmp = max(zs); % levee crest elevation relative to sea level
            if tmp_depth == 0 % channel depth is sometimes 0 at the entrance
                depthmax = [depthmax, NaN]; % if this is the case we want all variables to be NaN
                depth = [depth, NaN];
                zc = [zc;NaN];
                zlev = [zlev;NaN];
            else 
                depthmax = [depthmax, tmp_depth]; % otherwise we want to append data to depth array
                depth = [depth,tmp_depth]; % to save all channel depths for each radial transect
                zc = [zc;zc_tmp]; % to save channel bed elevations
                zlev = [zlev;zlev_tmp]; % to save channel levee elevations
            end
        end
        depthmax = depthmax'; % we don't want to average in the non-trunk channels, as they are often a lot shallower
        [val, idx] = max(depthmax); % find index of max depth
        % save channel depth for each distance through time
        if isempty(depthmax) % no depths?
            Hc_18(k,i) = NaN; % save NaN
            zc_18(k,i) = NaN;
            zlev_18(k,i) = NaN;
        else
            Hc_18(k,i) = val; % value of max depth
            zc_18(k,i) = zc(idx); % corresponding channel bed elevation
            zlev_18(k,i) = zlev(idx); % corresponding levee elevation
        end

        Hc_all18 = [Hc_all18;depth']; % compile all channel depths
        % now we will calculate aggradation rates
        % total aggradation
        dz_t = dz(~isnan(dz(:)));
        if isempty(dz_t)
            agg18_mean(k,i) = NaN;
            agg18_median(k,i) = NaN;
        else
            agg18_mean(k,i) = mean(dz_t/2, 'omitnan'); % mm/hr
            agg18_median(k,i) = median(dz_t/2, 'omitnan'); % mm/hr
        end
        % channel aggradation
        dz_c = dz_chan(~isnan(dz_chan(:)));
        if isempty(dz_c)
            c_agg18_mean(k,i) = NaN;
            c_agg18_median(k,i) = NaN;
        else 
            c_agg18_mean(k,i) = mean(dz_c/2, 'omitnan'); % mm/hr
            c_agg18_median(k,i) = median(dz_c/2, 'omitnan'); % mm/hr
        end
        % far-field aggradation
        dz_ff = dz_ff(~isnan(dz_ff(:)));
        if isempty(dz_ff)
            ff_agg18_mean(k,i) = NaN;
            ff_agg18_median(k,i) = NaN;
        else
            ff_agg18_mean(k,i) = mean(dz_ff/2, 'omitnan'); % mm/hr
            ff_agg18_median(k,i) = median(dz_ff/2, 'omitnan'); % mm/hr
        end
        % difference in agg between channel and ff
        if isempty(dz_c) || isempty(dz_ff)
            d_agg18_mean(k,i) = NaN;
            d_agg18_median(k,i) = NaN;
        else
            d_agg18_mean(k,i) = mean(dz_c/2, 'omitnan') - mean(dz_ff/2, 'omitnan'); % mm/hr
            d_agg18_median(k,i) = median(dz_c/2, 'omitnan') - median(dz_ff/2, 'omitnan'); % mm/hr
        end
    end
    % save depth percentiles for each radial segment
    for jj = 1:1:100
        if isempty(depth)
            Hc_tot18(k,jj) = NaN;
        else
            Hc_tot18(k,jj) = prctile(depth,jj);
        end
    end
end 

% treatment
% initialize empty matrices to fill 
rad_dist = 0:5:3100; % each radial transect is 50 mm or 5 cm
Hc_19 = NaN(length(rad_dist),nt_19-1); % channel depth with distance and time to use in aggradation below
Hc_tot19 = NaN(length(rad_dist),100); % to save percentile channel depth with distance
Hc_all19 = []; % all channel depths
zff_19 = NaN(length(rad_dist),nt_19-1); % to save mean far-field elevation rsl
zc_19 = NaN(length(rad_dist),nt_19-1); % to save channel bed elevation rsl
zlev_19 = NaN(length(rad_dist),nt_19-1); % to save levee crest elevation rsl
agg19_mean = NaN(length(rad_dist),nt_19-1); % total aggradation rate for distance through time
c_agg19_mean = NaN(length(rad_dist),nt_19-1); % channel aggradation rate for distance through time
ff_agg19_mean = NaN(length(rad_dist),nt_19-1); % far-field aggradation rate for distance through time
d_agg19_mean = NaN(length(rad_dist),nt_19-1); % difference in mean aggradation between channel and ff for distance through time
agg19_median = NaN(length(rad_dist),nt_19-1); % total aggradation rate for distance through time
c_agg19_median = NaN(length(rad_dist),nt_19-1); % channel aggradation rate for distance through time
ff_agg19_median = NaN(length(rad_dist),nt_19-1); % far-field aggradation rate for distance through time
d_agg19_median = NaN(length(rad_dist),nt_19-1); % difference in mean aggradation between channel and ff for distance through time
% loop through radial distances 
for k = 2:(length(rad_dist)-1) % loop to run through different radial distances from the end of the entrance channel.
    % print which radial segment the loop is on
    caption = sprintf('Radial segment %d', k);
    caption2 = sprintf(' of %d', length(rad_dist)-1);
    fprintf('%s\n', strcat(caption,caption2))

    % section to find x,y nodes for radial transect and generate matrix of
    % cross section topo and channel (yes/no) data
    idx = dd19(dd19 >= rad_dist(k) & dd19 < rad_dist(k+1));
    radial_dd = dd19 >= rad_dist(k) & dd19 < rad_dist(k+1); % can look at this using imagesc(radial_dd) to visualize a 0.1 m radial transect
    % initialize empty matrices 
    depth = [];
    for i = 1:(size(CM_buffer_19,3))
        % only want area on the delta top >-9 mm rsl or in the channels 
        CMtmp = CM_buffer_19(:,:,i); % we need tmp channel map
        CMtmp(isnan(CMtmp))=0; % to turn NaNs back to 0
        area = (terr_area19(:,:,i) + CMtmp).*basin19; %so we can add to area
        area(area>1)=1; % channel or >-9mm rsl = 1
        area(area==0)=NaN; % turn 0s back to NaN
        
        % dz
        rads = double(radial_dd); % turn into double
        rads(rads==0)=NaN; % so we can turn 0s to ones
        dz = dZD_19(:,:,i).*rads.*area; % so we can get dz only in the radius
        
        % channel and ff data in radial transect
        dz_chan = dZD_19_chan_buffer(:,:,i).*rads.*area; % multiply by radial transect 
        dz_ff = dZD_19_FF_buffer(:,:,i).*rads.*area; % multiply by radial transect

        % elevation data for channel depths
        z = zrsl19(:,:,i); % elevation data for channel depth
        is_shot = CM_buffer_19(:,:,i).*rads.*area; % in channel or no?
        z = z.*is_shot; % elevation in the channels
        is_shot(isnan(is_shot)) = 0; % we can't have NaNs in bwlabel
        
        % mean far-field elevation relative to sea level
        ztmp = zrsl19(:,:,i).*rads.*area.*FF_buffer_19(:,:,i); % far-field elevations
        if sum(~isnan(ztmp(:)))==0 % if no far-field in transect
            zff_19(k,i) = NaN; % save NaN
        else
            zff_19(k,i) = median(ztmp(:), 'omitnan'); % otherwise save median elevation of far-field transect
        end

        % how many individual channel segments ?
        [label,n] = bwlabel(is_shot); % label gives a unique number to each individual segement recognized; n is number of channels (or segments)
       
        % channel depth
        % we will calculate depth for each channel here
        depthmax = [];
        zc = [];
        zlev = [];
        for j = 1:n % loop through channel segments
            tmp = label;
            tmp(tmp ~= j) = NaN; % remove all data that is not in segement n
            tmp(tmp > 0) = 1; % turn segement data to 1
            tmp(tmp < 1) = NaN; % everything else is NaN 
            z_tmp = z.*tmp; % segment elevations to get depth
            zs = z_tmp(~isnan(z_tmp)); % remove data not in the channel
            tmp_depth = (max(zs) - min(zs)); % channel depth in mm; we will take the max elevation as levee elevation and min elevation as thalweg                  
            zc_tmp = min(zs); % channel bed elevation relative to sea level
            zlev_tmp = max(zs); % levee crest elevation relative to sea level
            if tmp_depth == 0 % sometimes channel depth is 0 near the entrance
                depthmax = [depthmax, NaN]; % so append NaN
                depth = [depth, NaN];
                zc = [zc;NaN];
                zlev = [zlev;NaN];
            else 
                depthmax = [depthmax, tmp_depth]; % otherwise append depths
                depth = [depth,tmp_depth]; % to save all channel depths for each radial transect
                zc = [zc;zc_tmp]; % channel bed elevation
                zlev = [zlev;zlev_tmp]; % levee elevation
            end
        end
        depthmax = depthmax'; % we don't want to average in the non-trunk channels, as they are often a lot shallower
        [val, idx] = max(depthmax); % find index of deepest channel
        % save channel depth for each distance through time
        if isempty(depthmax) % if no channels in transect
            Hc_19(k,i) = NaN;
            zc_19(k,i) = NaN;
            zlev_19(k,i) = NaN;
        else
            Hc_19(k,i) = val; % otherwise save value of max channel depth
            zc_19(k,i) = zc(idx); % and corresponding channel bed elevation
            zlev_19(k,i) = zlev(idx); % and corresponding levee elevation
        end

        Hc_all19 = [Hc_all19;depth']; % compile all channel depths
        % now we will calculate aggradation rates
        % total aggradation
        dz_t = dz(~isnan(dz(:)));
        if isempty(dz_t)
            agg19_mean(k,i) = NaN;
            agg19_median(k,i) = NaN;
        else
            agg19_mean(k,i) = mean(dz_t/2, 'omitnan'); % mm/hr
            agg19_median(k,i) = median(dz_t/2, 'omitnan'); % mm/hr
        end
        % channel aggradation
        dz_c = dz_chan(~isnan(dz_chan(:)));
        if isempty(dz_c)
            c_agg19_mean(k,i) = NaN;
            c_agg19_median(k,i) = NaN;
        else 
            c_agg19_mean(k,i) = mean(dz_c/2, 'omitnan'); % mm/hr
            c_agg19_median(k,i) = median(dz_c/2, 'omitnan'); % mm/hr
        end
        % far-field aggradation
        dz_ff = dz_ff(~isnan(dz_ff(:)));
        if isempty(dz_ff)
            ff_agg19_mean(k,i) = NaN;
            ff_agg19_median(k,i) = NaN;
        else
            ff_agg19_mean(k,i) = mean(dz_ff/2, 'omitnan'); % mm/hr
            ff_agg19_median(k,i) = median(dz_ff/2, 'omitnan'); % mm/hr
        end
        % difference in agg between channel and ff
        if isempty(dz_c) || isempty(dz_ff)
            d_agg19_mean(k,i) = NaN;
            d_agg19_median(k,i) = NaN;
        else
            d_agg19_mean(k,i) = mean(dz_c/2, 'omitnan') - mean(dz_ff/2, 'omitnan'); % mm/hr
            d_agg19_median(k,i) = median(dz_c/2, 'omitnan') - median(dz_ff/2, 'omitnan'); % mm/hr
        end
    end
    % save depth percentiles for each radial segment
    for jj = 1:1:100
        if isempty(depth)
            Hc_tot19(k,jj) = NaN;
        else
            Hc_tot19(k,jj) = prctile(depth,jj);
        end
    end
end

%% Remove outliers from channel depth data if you want
mean(Hc_all18(:), 'omitnan')
std(Hc_all18(:), 'omitnan')
mean(Hc_all19(:), 'omitnan')
std(Hc_all19(:), 'omitnan')
median(Hc_all18(:), 'omitnan')
median(Hc_all19(:), 'omitnan')

%% Calculate compensation timescale and channel in-filling rate
hc = Hc_18(:,1:279);
hc18 = mean(hc,2, 'omitnan');
agg18_tmp = mean(agg18_mean,2, 'omitnan');
dagg18 = mean(d_agg18_mean,2, 'omitnan');
Tc_18 = hc18./agg18_tmp;
Ta_18 = hc18./dagg18;

hc = Hc_19(:,1:279);
hc19 = mean(hc,2, 'omitnan');
agg19_tmp = mean(agg19_mean,2, 'omitnan');
dagg19 = mean(d_agg19_mean,2, 'omitnan');
Tc_19 = hc19./agg19_tmp;
Ta_19 = hc19./dagg19;
%infill_std = infill*(sqrt((dagg_std/dagg)^2+(depth_std18(k)/depth18(k))^2)); % error propagation
% control
% dagg_std = sqrt((std(c_agg18)^2)+(std(f_agg18)^2)/2); % error propagation

%% Plot channel depth
rad_dist = rad_dist./1000.;

% make figure
fig = figure();
plot(rad_dist, hc18', 'b-', 'LineWidth', 2)
hold on
plot(rad_dist, hc19', 'g-', 'LineWidth', 2)
ylabel('trunk channel depth (mm)')
xlabel('distance from apex (m)') 
ylim([0 35])
legend('control', 'treatment')%, 'control IQR', 'treatment IQR')
set(gca, 'YMinorTick', 'On', 'XMinorTick', 'On', 'XAxisLocation', 'bottom', 'XAxisLocation', 'bottom')
set(gcf, 'PaperUnits', 'inches');
y_width=7.25;x_width=9.125;
set(gcf, 'PaperPosition', [0 0 x_width y_width]);
saveas(fig, '../figures/esurf_Figure6b.pdf')

%% Plot the data: compensation timescale
% create arrays
tc18 = mean(Tc_18,2, 'omitnan');
chan_array18 = [rad_dist; tc18'];
cols = any(isnan(chan_array18),1);
chan_array18(:,cols) = [];

tc19 = mean(Tc_19,2, 'omitnan');
chan_array19 = [rad_dist; tc19'];
cols = any(isnan(chan_array19),1);
chan_array19(:,cols) = [];

y18 = chan_array18(2,:); % your mean vector;
x18 = chan_array18(1,:);

y19 = chan_array19(2,:); % your mean vector;
x19 = chan_array19(1,:);

% make figure
fig = figure();
plot(x18, y18, 'b-', 'LineWidth', 2)
hold on
plot(x19, y19, 'g-', 'LineWidth', 2)
ylabel('compensation timescale (hrs)')
ylim([0 120])
legend('control mean', 'treatment mean')
set(gca, 'XMinorTick', 'On', 'YMinorTick', 'On')
set(gcf, 'PaperUnits', 'inches');
y_width=7.25;x_width=9.125;
set(gcf, 'PaperPosition', [0 0 x_width y_width]);
saveas(fig, '../figures/esurf_Figure6X.pdf')

%% Plot the data: levee elevation, channel elevation, far-field elevation
% create arrays
% channel
zchan18 = mean(zc_18,2, 'omitnan');
stdev18 = std(zc_18,[],2, 'omitnan');
chan_array18 = [rad_dist; zchan18'; stdev18'];
cols = any(isnan(chan_array18),1);
chan_array18(:,cols) = [];

zchan19 = mean(zc_19,2, 'omitnan');
stdev19 = std(zc_19,[],2, 'omitnan');
chan_array19 = [rad_dist; zchan19'; stdev19'];
cols = any(isnan(chan_array19),1);
chan_array19(:,cols) = [];

%fill standard deviation
yc18 = chan_array18(2,:); % your mean vector;
xc18 = chan_array18(1,:);
std18 = chan_array18(3,:);
curve1c_18 = yc18 + std18;
curve2c_18 = yc18 - std18;

yc19 = chan_array19(2,:); % your mean vector;
xc19 = chan_array19(1,:);
std19 = chan_array19(3,:);
curve1c_19 = yc19 + std19;
curve2c_19 = yc19 - std19;

% create arrays
% levee
zl18 = mean(zlev_18,2, 'omitnan');
stdev18 = std(zlev_18,[],2, 'omitnan');
chan_array18 = [rad_dist; zl18'; stdev18'];
cols = any(isnan(chan_array18),1);
chan_array18(:,cols) = [];

zl19 = mean(zlev_19,2, 'omitnan');
stdev19 = std(zlev_19,[],2, 'omitnan');
chan_array19 = [rad_dist; zl19'; stdev19'];
cols = any(isnan(chan_array19),1);
chan_array19(:,cols) = [];

%fill standard deviation
yl18 = chan_array18(2,:); % your mean vector;
xl18 = chan_array18(1,:);
std18 = chan_array18(3,:);
curve1l_18 = yl18 + std18;
curve2l_18 = yl18 - std18;

yl19 = chan_array19(2,:); % your mean vector;
xl19 = chan_array19(1,:);
std19 = chan_array19(3,:);
curve1l_19 = yl19 + std19;
curve2l_19 = yl19 - std19;

% create arrays
%far field 
zf18 = mean(zff_18,2, 'omitnan');
stdev18 = std(zff_18,[],2, 'omitnan');
chan_array18 = [rad_dist; zf18'; stdev18'];
cols = any(isnan(chan_array18),1);
chan_array18(:,cols) = [];

zf19 = mean(zff_19,2, 'omitnan');
stdev19 = std(zff_19,[],2, 'omitnan');
chan_array19 = [rad_dist; zf19'; stdev19'];
cols = any(isnan(chan_array19),1);
chan_array19(:,cols) = [];

%fill standard deviation
yf18 = chan_array18(2,:); % your mean vector;
xf18 = chan_array18(1,:);
std18 = chan_array18(3,:);
curve1f_18 = yf18 + std18;
curve2f_18 = yf18 - std18;

yf19 = chan_array19(2,:); % your mean vector;
xf19 = chan_array19(1,:);
std19 = chan_array19(3,:);
curve1f_19 = yf19 + std19;
curve2f_19 = yf19 - std19;

% make figure
fig = figure();
plot(xc18, yc18, 'b-', 'LineWidth', 2)
hold on
plot(xc19, yc19, 'g-', 'LineWidth', 2)
plot(xl18, yl18, 'b--', 'LineWidth', 2)
plot(xl19, yl19, 'g--', 'LineWidth', 2)
plot(xf18, yf18, 'b:', 'LineWidth', 2)
plot(xf19, yf19, 'g:', 'LineWidth', 2)
ylabel('elevation relative to sea level (mm)')
legend('control channel', 'treatment channel', 'control levee', 'treatment levee', 'control far-field', 'treatment far-field')
set(gca, 'XMinorTick', 'On', 'YMinorTick', 'On')
set(gcf, 'PaperUnits', 'inches');
y_width=7.25;x_width=9.125;
set(gcf, 'PaperPosition', [0 0 x_width y_width]);
saveas(fig, '../figures/esurf_Figure6X.pdf')


%% Calculate superelevatin
Lr18 = mean(zlev_18,2,'omitnan')-mean(zff_18,2,'omitnan');
Lr19 = mean(zlev_19,2,'omitnan')-mean(zff_19,2,'omitnan');
h18 = mean(Hc_18,2,'omitnan');
h19 = mean(Hc_19,2,'omitnan');

SER18 = Lr18./h18;
SER19 = Lr19./h19;

fig = figure;
plot(rad_dist, SER18', 'b-', 'LineWidth', 2)
hold on
plot(rad_dist, SER19', 'g-', 'LineWidth', 2)
ylim([0 10])
ylabel('superelevation ratio (-)')
legend('control', 'treatment')
set(gca, 'XMinorTick', 'On', 'YMinorTick', 'On')
%% Plot the data: avulsion timescale
% create arrays
ta18 = mean(Ta_18,2, 'omitnan');
chan_array18 = [rad_dist; ta18'];
cols = any(isnan(chan_array18),1);
chan_array18(:,cols) = [];

ta19 = mean(Ta_19,2, 'omitnan');
chan_array19 = [rad_dist; ta19'];
cols = any(isnan(chan_array19),1);
chan_array19(:,cols) = [];

%clean data near apex
idx = chan_array18(1, :) < 0.1;
chan_array18(:, idx) = NaN;

idx = chan_array19(1, :) < 0.1;
chan_array19(:, idx) = NaN;

%fill standard deviation
y18 = chan_array18(2,:); % your mean vector;
x18 = chan_array18(1,:);

y19 = chan_array19(2,:); % your mean vector;
x19 = chan_array19(1,:);

fig = figure();
plot(x18, y18, 'b-', 'LineWidth', 2)
hold on
plot(x19, y19, 'g-', 'LineWidth', 2)
ylabel('avulsion timescale (hrs)')
xlabel('distance from apex (m)')
ylim([0 150])
legend('control mean', 'treatment mean')
set(gca, 'XMinorTick', 'On', 'YMinorTick', 'On')
set(gcf, 'PaperUnits', 'inches');
y_width=7.25;x_width=9.125;
set(gcf, 'PaperPosition', [0 0 x_width y_width]);
saveas(fig, '../figures/esurf_Figure6d.pdf')
%% Plot the data: Figure 5
% error bars for plotting
% channel aggradation
chan_agg_mean_18 = mean(c_agg18_mean,2, 'omitnan');
chan_agg_std_18 = mean(c_agg18_mean,2, 'omitnan');
chan_array18 = [rad_dist; chan_agg_mean_18'; chan_agg_std_18'];
cols = any(isnan(chan_array18),1);
chan_array18(:,cols) = [];

chan_agg_mean_19 = mean(c_agg19_mean,2, 'omitnan');
chan_agg_std_19 = mean(c_agg19_mean,2, 'omitnan');
chan_array19 = [rad_dist; chan_agg_mean_19'; chan_agg_std_19'];
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
ff_agg_mean_18 = mean(ff_agg18_mean,2, 'omitnan');
ff_agg_std_18 = mean(ff_agg18_mean,2, 'omitnan');
ff_array18 = [rad_dist; ff_agg_mean_18'; ff_agg_std_18'];
cols = any(isnan(ff_array18),1);
ff_array18(:,cols) = [];

ff_agg_mean_19 = mean(ff_agg19_mean,2, 'omitnan');
ff_agg_std_19 = mean(ff_agg19_mean,2, 'omitnan');
ff_array19 = [rad_dist; ff_agg_mean_19'; ff_agg_std_19'];
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
d_agg_mean_18 = mean(d_agg18_mean,2, 'omitnan');
d_agg_std_18 = mean(d_agg18_mean,2, 'omitnan');
dd_array18 = [rad_dist; d_agg_mean_18'; d_agg_std_18'];
cols = any(isnan(dd_array18),1);
dd_array18(:,cols) = [];

d_agg_mean_19 = mean(d_agg19_median,2, 'omitnan');
d_agg_std_19 = mean(d_agg19_median,2, 'omitnan');
dd_array19 = [rad_dist; d_agg_mean_19'; d_agg_std_19'];
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
ylim([0 3])
xlim([0 3])
ylabel('aggradation rate (mm/hr)')
xlabel('distance from apex (m)')
legend('control channel', 'treatment channel', 'control far-field', 'treatment far-field', 'RSLR{_b}')
set(gca,'XMinorTick','on','YMinorTick','on')
y_width=7.25;x_width=9.125;
set(gcf, 'PaperPosition', [0 0 x_width y_width]);
saveas(fig, '../figures/esurf_Figure6a.pdf')

% Figure 5b
fig = figure();
plot(dd_x18, dd_y18, 'b-', 'LineWidth', 2)
hold on
plot(dd_x19, dd_y19, 'g-', 'LineWidth', 2)
ylabel('channel in-filling rate (mm/hr)')
ylim([0 3])
xlabel('distance from apex (m)') 
legend('control mean', 'treatment mean') 
set(gca,'XMinorTick','on','YMinorTick','on')
y_width=7.25;x_width=9.125;
set(gcf, 'PaperPosition', [0 0 x_width y_width]);
saveas(fig, '../figures/esurf_Figure6c.pdf')
