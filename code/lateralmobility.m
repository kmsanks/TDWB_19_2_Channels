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

%% Find fraction of delta that is above sea level 
%control
pland_18 = [];
t_18 = [];
for i = 1:nt_18
    t_18 = [t_18;i*dt_18]; %time (hours)
    z = ZD_18(:,:,i); %elevation (mm)
    z = z - (baselevel_rr*i*dt_18+ocean_zero); %elevation relative to sea level
    z(z > 0) = 1; %above sea level binary of 1s
    z(z < 0) = 0; %below sea level binary of 0s
    pland_18(:,:,i) = z; %save land mask
end
%treatment
pland_19 = [];
t_19 = [];
for i = 1:nt_19
    t_19 = [t_19;(i-1)*dt_19]; %time (hours)
    z = ZD_19(:,:,i); %elevation (mm)
    z = z - (baselevel_rr*(i-1)*dt_19+ocean_zero); %elevation relative to sea level
    z(z > 0) = 1; %above sea level binary of 1s
    z(z < 0) = 0; %below sea level binary of 0s
    pland_19(:,:,i) = z; %save land mask
end 

%control
pland_18_sum = sum(pland_18,3); %how many times each grid cell is above sea level
pland_18_sum(pland_18_sum == 0) = NaN; %not on land make NaN
pland_18_frac = pland_18_sum/560; %fraction of time pixel is on land

frac_18 = [];
area_pland_18 = []; %area on land for various amounts of the experiment
for i = round(linspace(0,1,101),2)
    frac_18 = [frac_18;i];
    area = pland_18_frac(pland_18_frac >= i); %& pland_18_frac < i+1
    area(area > 0) = 1;
    areasum = sum(area);
    area_pland_18_i = (areasum)*2.5*10^-5; %area of 1 pixel
    area_pland_18 = [area_pland_18; area_pland_18_i];     
end 

%treatment
pland_19_sum = sum(pland_19,3);
pland_19_sum(pland_19_sum == 0) = NaN;
pland_19_frac = pland_19_sum/281;

frac_19 = [];
area_pland_19 = [];
for i = round(linspace(0,1,101),2)
    frac_19 = [frac_19;i];
    area = pland_19_frac(pland_19_frac >= i); %& pland_18_frac < i+1
    area(area > 0) = 1;
    areasum = sum(area);
    area_pland_19_i = (areasum)*2.5*10^-5; %area of 1 pixel
    area_pland_19 = [area_pland_19; area_pland_19_i];     
end 

% Find the 50% area which will be used to calculate lateral mobility
% statistics
frac_area_array = [frac_18, area_pland_18];
frac50_index = find(frac_area_array(:,1) == 0.50); %50 percent area
area_frac50_18 = frac_area_array(frac50_index,2);

frac_area_array_19 = [frac_19, area_pland_19];
frac50_index_19 = find(frac_area_array_19(:,1) == 0.50); %50 percent area
area_frac50_19 = frac_area_array_19(frac50_index_19,2);

%% Plot the area above sea level for different amount of time: SI Fig. B1
fig = figure();
plot(frac_18, area_pland_18, 'bo');
hold on
plot(frac_19, area_pland_19, 'go');
xline(0.5, 'LineWidth',2);
yline(area_frac50_18, 'b-', 'LineWidth',2);
yline(area_frac50_19, 'g-', 'LineWidth',2);
xlabel('f_{land} (-)');
ylabel('area (m^{2})');
legend('control', 'treatment', '50%', 'control 50%', 'treatment 50%');
set(gca, 'XMinorTick','on', 'YMinorTick','On')
y_width=7.25;x_width=9.125
set(gcf, 'PaperPosition', [0 0 x_width y_width]);
saveas(fig, '../figures/esurf_FigureB1.pdf')

%% Now we will calculate the lateral channel mobility
% Here will will determine how long it takes for 90% of the delta to be
% visited by a channel through time 

% control
% crop channels to 50% land area
pland_18_frac_nan = pland_18_frac;
pland_18_frac_nan(pland_18_frac_nan < 0.5) = NaN;
pland_18_frac_nan(pland_18_frac_nan > 0) = 1;
CM_18_crop = [];
for i = 1:((size(CM_18,3)))
    CM_18_nan = CM_18(:,:,i).*pland_18_frac_nan; % crop channel maps to 50% land area
    CM_18_crop(:,:,i) = CM_18_nan; % save new channel maps
end

% calculate lateral channel mobility
area = area_frac50_18; %0.5 land above sea level 2.1868 m^2
lat_chan_mob_18 = [];
for i = 1:((size(CM_18_crop,3)-60)) % loop through time 
    for ii = i:size(CM_18_crop,3) % loop through time
        if ii == i % if same time step, then area covered by channel is just current area covered
            CM_n = CM_18_crop(:,:,ii); % channel map at time ii
            lat_chan_mob_18_i = 1-((sum(CM_n(:), 'omitnan')*2.5*10^-5)/area); %area covered is = channel area
            hour_18_i = 0;
        else
            CM_n = CM_n + CM_18_crop(:,:,ii); %add all channel maps together from time i to ii
            CM_n(CM_n > 1) = 1; %if it has ever been a channel, turn to 1
            lat_chan_mob_18_i = 1-((sum(CM_n(:), 'omitnan')*2.5*10^-5)/area); %area covered = all channel area from time i to ii
            hour_18_i = ii-i; % how much time has passed between i and ii?
            lat_chan_mob_18_i = [hour_18_i, lat_chan_mob_18_i];
            lat_chan_mob_18 = [lat_chan_mob_18;lat_chan_mob_18_i];
        end
    end
end 

% treatment
% Replace timesteps with no channel maps with the channel map from the
% next time step for the treatment experiment
for i = (560-1):-1:1 %I know 560 has a channel map, so I can start at 279 and replace with channel map that comes after, this will work for ones that have multiple no maps in a row
    if nansum(nansum(CM_19(:,:,i))) == 0
        CM_19(:,:,i) = CM_19(:,:,i+1);
    end
end
% now crop the channel maps to the 50% land area for the treatment
pland_19_frac_nan = pland_19_frac;
pland_19_frac_nan(pland_19_frac_nan < 0.5) = NaN;
pland_19_frac_nan(pland_19_frac_nan > 0) = 1;
CM_19_crop = [];
for i = 1:((size(CM_19,3)))
    CM_19_nan = CM_19(:,:,i).*pland_19_frac_nan;
    CM_19_crop(:,:,i) = CM_19_nan;
end

area = area_frac50_19; %0.5 land above sea level 1.9329 m^2
lat_chan_mob_19 = [];
for i = 1:((size(CM_19_crop,3)-60))
    for ii = i:size(CM_19_crop,3)
        if ii == i
            CM_n = CM_19_crop(:,:,ii);
            lat_chan_mob_19_i = 1-((sum(CM_n(:), 'omitnan')*2.5*10^-5)/area);
            hour_19_i = 0;
        else
            CM_n = CM_n + CM_19_crop(:,:,ii);
            CM_n(CM_n > 1) = 1;
            lat_chan_mob_19_i = 1-((sum(CM_n(:), 'omitnan')*2.5*10^-5)/area);
            hour_19_i = ii-i;
            lat_chan_mob_19_i = [hour_19_i, lat_chan_mob_19_i];
            lat_chan_mob_19 = [lat_chan_mob_19;lat_chan_mob_19_i];
        end
    end
end 

%% Calculate lateral mobility through time for plotting and statistics
% control
[ud,~,iy] = unique(lat_chan_mob_18(:,1));  
mean_lat_chan_mob_18 = [ud, accumarray(iy,lat_chan_mob_18(:,2),[],@mean)];
std_lat_chan_mob_18 = [ud, accumarray(iy,lat_chan_mob_18(:,2),[],@std)];
% treatment
[ud,~,iy] = unique(lat_chan_mob_19(:,1));  
mean_lat_chan_mob_19 = [ud, accumarray(iy,lat_chan_mob_19(:,2),[],@mean)];
std_lat_chan_mob_19 = [ud, accumarray(iy,lat_chan_mob_19(:,2),[],@std)];
std_lat_chan_mob_19(std_lat_chan_mob_19 < 0) = 0;

%differnce to plot contours on each other
xdif = xentrance_19 - xentrance_18;
ydif = yentrance_19 - yentrance_18;

%subtract xdif and ydif from 19 matrix
pland_19_frac_shift = pland_19_frac(xdif:end, ydif:end, :);
pland_19_frac_nan_shift = pland_19_frac_nan(xdif:end, ydif:end, :);

% plot contours for control and treatment pland50
% v = [0.5,0.5];
% figure()
% imagesc(pland_18_frac_nan)
% contour(pland_18_frac,v, 'b-')
% hold on
% %imagesc(pland_19_frac_nan_shift)
% %alpha(pland_19_frac_nan_shift, 'clear', 0.1)
% contour(pland_19_frac_shift, v, 'g-')
% legend('control', 'treatment')
% set(gca, 'Ydir', 'reverse');

%% Now we will calculate the fraction unmodified statistic (f_um)
% or how long it takes the delta to accumulate 1mm of sediment

% control
% crop elevations to 50% area mask
% Here we will only use every other hour for the control experiment to
% avoid Saddler effects
ZD_18_crop = ZD_18(:,:,1:2:end).*pland_18_frac_nan;
ZD_19_crop = ZD_19.*pland_19_frac_nan;

nt=280;
area = area_frac50_18; %0.5 land above sea level 1.9329 m^2
f_um_18 = [];
for i = 1:nt-1
    for ii = i:nt
        dz_18_i = ZD_18_crop(:,:,ii)-ZD_18_crop(:,:,i);
        dz_18_i(dz_18_i < 1) = 0;
        dz_18_i(dz_18_i > 1) = 1;
        f_um_18_i = 1-(sum(sum(dz_18_i(:), 'omitnan')*2.5*10^-5)/area);
        ts_i = (ii-i)*2; % hours
        f_um_18_i = [ts_i, f_um_18_i];
        f_um_18 = [f_um_18;f_um_18_i];
    end
end

nt=281;
area = area_frac50_19; %0.5 land above sea level 1.9329 m^2
f_um_19 = [];
for i = 1:(size(ZD_19_crop,3)-1)
    for ii = i:size(ZD_19_crop,3)
        dz_19_i = ZD_19_crop(:,:,ii)-ZD_19_crop(:,:,i);
        dz_19_i(dz_19_i < 1) = 0;
        dz_19_i(dz_19_i > 1) = 1;
        f_um_19_i = 1-(sum(sum(dz_19_i(:), 'omitnan')*2.5*10^-5)/area);
        ts_i = (ii-i)*2; % hours
        f_um_19_i = [ts_i, f_um_19_i];
        f_um_19 = [f_um_19;f_um_19_i];
    end
end

[ud,~,iy] = unique(f_um_18(:,1));  
mean_f_um_18 = [ud, accumarray(iy,f_um_18(:,2),[],@mean)];
std_f_um_18 = [ud, accumarray(iy,f_um_18(:,2),[],@std)];

[ud,~,iy] = unique(f_um_19(:,1));  
mean_f_um_19 = [ud, accumarray(iy,f_um_19(:,2),[],@mean)];
std_f_um_19 = [ud, accumarray(iy,f_um_19(:,2),[],@std)];

% plot with standard deviation
% lateral channel mobility
array18_lat = [(1:559); mean_lat_chan_mob_18(:,2)'; std_lat_chan_mob_18(:,2)'];
cols = any(isnan(array18_lat),1);
array18_lat(:,cols) = [];

array19_lat = [(1:559);mean_lat_chan_mob_19(:,2)'; std_lat_chan_mob_19(:,2)'];
cols = any(isnan(array19_lat),1);
array19_lat(:,cols) = [];

% fill standard deviation
y18_lat = array18_lat(2,:); % your mean vector;
x18_lat = array18_lat(1,:);
std18_lat = array18_lat(3,:);
curve1_18_lat = y18_lat + std18_lat;
curve2_18_lat = y18_lat - std18_lat;

y19_lat = array19_lat(2,:); % your mean vector;
x19_lat = array19_lat(1,:);
std19_lat = array19_lat(3,:);
curve1_19_lat = y19_lat + std19_lat;
curve2_19_lat = y19_lat - std19_lat;

% fraction unmodified
array18_fum = [(1:2:560); mean_f_um_18(:,2)'; std_f_um_18(:,2)'];
cols = any(isnan(array18_fum),1);
array18_fum(:,cols) = [];

array19_fum = [(1:2:562);mean_f_um_19(:,2)'; std_f_um_19(:,2)'];
cols = any(isnan(array19_fum),1);
array19_fum(:,cols) = [];

% fill standard deviation
y18_fum = array18_fum(2,:); % your mean vector;
x18_fum = array18_fum(1,:);
std18_fum = array18_fum(3,:);
curve1_18_fum = y18_fum + std18_fum;
curve2_18_fum = y18_fum - std18_fum;

y19_fum = array19_fum(2,:); % your mean vector;
x19_fum = array19_fum(1,:);
std19_fum = array19_fum(3,:);
curve1_19_fum = y19_fum + std19_fum;
curve2_19_fum = y19_fum - std19_fum;

%% Plot the data: Figure 7ab
fig = figure();
%plot lateral channel mobility
subplot(1,2,1)
plot(x18_lat, y18_lat, 'b', 'LineWidth', 2)
hold on
plot(x19_lat, y19_lat, 'g', 'LineWidth', 2)
patch([x18_lat fliplr(x18_lat)], [curve1_18_lat fliplr(curve2_18_lat)], 'b')
patch([x19_lat fliplr(x19_lat)], [curve1_19_lat fliplr(curve2_19_lat)], 'g')
alpha(0.15)
plot(x18_lat, y18_lat, 'b', 'LineWidth', 2)
plot(x19_lat, y19_lat, 'g', 'LineWidth', 2)
ylim([0 1])
xlim([0 300])
axis('square')
grid on
grid minor
ylabel('fraction not visited by channel (-)')
xlabel('measurement window (hr)') 
legend('control mean', 'treatment mean', 'control stdev', 'treatment stdev')
%plot fraction unmodified 
subplot(1,2,2)
plot(x18_fum, y18_fum, 'b', 'LineWidth', 2)
hold on
plot(x19_fum, y19_fum, 'g', 'LineWidth', 2)
patch([x18_fum fliplr(x18_fum)], [curve1_18_fum fliplr(curve2_18_fum)], 'b')
patch([x19_fum fliplr(x19_fum)], [curve1_19_fum fliplr(curve2_19_fum)], 'g')
alpha(0.15)
plot(x18_fum, y18_fum, 'b', 'LineWidth', 2)
plot(x19_fum, y19_fum, 'g', 'LineWidth', 2)
ylim([0 1])
xlim([0 200])
grid on
grid minor
ylabel('fraction unmodified (-)')
xlabel('measurment window (hr)') 
legend('control mean', 'treatment mean', 'control stdev', 'treatment stdev')
axis('square')
set(gcf, 'PaperUnits', 'inches');
y_width=6;x_width=6;
set(gcf, 'PaperPosition', [0 0 x_width y_width]);
saveas(fig, '../figures/esurf_Figure7ab.pdf')

%% Plot the data: Figure B2
fig = figure();
subplot(1,2,1)
semilogy(mean_lat_chan_mob_18(:,1), mean_lat_chan_mob_18(:,2), 'bo');
hold on
semilogy(mean_lat_chan_mob_19(:,1), mean_lat_chan_mob_19(:,2), 'go');
grid on
grid minor
xlim([0 560]);
xlabel('measurement window (hr)');
ylabel('log(fraction not visited by channel (-))');
legend('control', 'treatment');
axis('square');
subplot(1,2,2)
semilogy(mean_f_um_18(:,1), mean_f_um_18(:,2), 'bo');
grid on
grid minor
hold on
semilogy(mean_f_um_19(:,1), mean_f_um_19(:,2), 'go');
grid on
grid minor
xlim([0 200]);
xlabel('measurement window (hr)');
ylabel('log(fraction unmodified (-))');
legend('control', 'treatment');
axis('square');
%patch([min(xlim) max(xlim) max(xlim) min(xlim)], [ybars(1) ybars(1), ybars(2) ybars(2)], [0.8 0.8 0.8])
set(gcf, 'PaperUnits', 'inches');
y_width=6;x_width=6;
set(gcf, 'PaperPosition', [0 0 x_width y_width]);
saveas(fig, '../figures/esurf_FigureB2.pdf')

%% Calculate timescales for lateral mobility and modification
% exponential fit for modification ts
x = mean_f_um_18(:,1);
x1 = x(1:60);
y = mean_f_um_18(:,2);
y1 = y(1:60);
f0 = ezfit(x1, y1, 'exp');
efold18 = -1/f0.m(2);

% exponential fit for modification ts
x = mean_f_um_19(:,1);
x2 = x(1:80);
y = mean_f_um_19(:,2);
y2 = y(1:80);
f1 = ezfit(x2, y2, 'exp');
efold19 = -1/f1.m(2);

% exponential fit for lateral channel mobility
x = mean_lat_chan_mob_18(:,1);
x1 = x(1:380);
y = mean_lat_chan_mob_18(:,2);
y1 = y(1:380);
f0 = ezfit(x1, y1, 'exp');
latfold18 = -1/f0.m(2);

% exponential fit for lateral mobility
x = mean_lat_chan_mob_19(:,1);
x2 = x(1:240);
y = mean_lat_chan_mob_19(:,2);
y2 = y(1:240);
f1 = ezfit(x2, y2, 'exp');
latfold19 = -1/f1.m(2);

%% Channel density and histogram
% How many hours does a channel occupy each delta pixel?
chan_dens18 = sum(CM_18_crop, 3, 'omitnan');
chan_dens19 = sum(CM_19_crop, 3, 'omitnan');

% Plot the data: Figure B6
fig = figure;
subplot(2,2,1)
imagesc(chan_dens18(1:600,109:509))
colorbar
caxis([0 250])
axis square
subplot(2,2,2)
imagesc(chan_dens19(120:719,214:614))
colorbar
caxis([0 250])
axis square
subplot(2,2,3)
histogram(chan_dens18(chan_dens18>0), 'Normalization', 'pdf', 'BinWidth',5, 'FaceColor', 'blue')
xline(median(chan_dens18(chan_dens18>0)), 'k-', 'linewidth', 2)
legend('control','median')
xlim([0 250])
ylim([0 .03])
subplot(2,2,4)
histogram(chan_dens19(chan_dens19>0), 'Normalization', 'pdf', 'BinWidth',5, 'FaceColor', 'green')
xline(median(chan_dens19(chan_dens19>0)), 'k-', 'linewidth', 2)
legend('treatment','median')
xlim([0 250])
ylim([0 .03])
y_width=7.25;x_width=9.125;
set(gcf, 'PaperPosition', [0 0 x_width y_width]);
saveas(fig, '../figures/esurf_FigureB6.pdf')

%% Consecutive channelization
% Now we will determine the longest consecutive time that each pixel remains a channel
% Note that 19 has some missing channel maps, so we want to skip these and
% not restart our counter here
% consec18 = nandiff(chanMaps_18_crop, 3);
% chan_dens19 = nansum(chanMaps_19_crop, 3);

max_consec_channel_18 = NaN(size(CM_18_nan));
for i = 1:size(CM_18,1)
    for j = 1:size(CM_18,2)
        consec = CM_18(i,j,:);
        consec = consec(:);
        consec = consec'; 
        n = find(diff(consec));
        n_consec = [n numel(consec)] - [0 n];
        bin = consec(n);
        bin(size(bin,2)+1) = consec(560);
        consec_array = cat(1,n_consec,bin);
        idx = find(consec_array(2,:)==1);
        if isempty(idx)
            max = NaN;
        else
            max = nanmax(consec_array(1,idx));
        end
        max_consec_channel_18(i,j) = max;
    end
end 

max_consec_channel_19 = NaN(size(CM_19_nan));
for i = 1:size(CM_19,1)
    for j = 1:size(CM_19,2)
        consec = CM_19(i,j,:);
        consec = consec(:);
        consec = consec'; 
        n = find(diff(consec));
        n_consec = [n numel(consec)] - [0 n];
        bin = consec(n);
        bin(size(bin,2)+1) = consec(560);
        consec_array = cat(1,n_consec,bin);
        idx = find(consec_array(2,:)==1);
        if isempty(idx)
            max = NaN;
        else
            max = nanmax(consec_array(1,idx));
        end
        max_consec_channel_19(i,j) = max;
    end
end

% now we will normalize by the total channelized time
chan_dens18 = sum(CM_18, 3, 'omitnan');
chan_dens19 = sum(CM_19, 3, 'omitnan');

max_dens18 = max_consec_channel_18./chan_dens18;
max_dens19 = max_consec_channel_19./chan_dens19;

% Plot the data: Figure B7
fig = figure();
subplot(3,2,1)
imagesc(max_consec_channel_18(1:600,109:509))
colorbar
caxis([0 100])
axis square
title('maximum consecutive channelized time (hrs; control)')
subplot(3,2,2)
imagesc(max_consec_channel_19(120:719,214:614))
colorbar
caxis([0 100])
axis square
title('maximum consecutive channelized time (hrs; treatment)')
subplot(3,2,3)
imagesc(max_dens18(1:600,109:509))
colorbar
caxis([0 1.])
axis square
title('channelizied time (-; control)')
subplot(3,2,4)
imagesc(max_dens19(120:719,214:614))
colorbar
caxis([0 1.])
axis square
title('channelizied time (-; treatment)')
subplot(3,2,5)
histogram(max_consec_channel_18(max_consec_channel_18 > 0), 'Normalization', 'pdf', 'BinWidth',2, 'FaceColor', 'blue')
xline(median(max_consec_channel_18(max_consec_channel_18>0)), 'k-', 'linewidth', 2)
legend('control','median')
xlabel('max consecutive hours channelized')
xlim([0 50])
ylim([0 0.15])
subplot(3,2,6)
histogram(max_consec_channel_19(max_consec_channel_19 > 0), 'Normalization', 'pdf', 'BinWidth',2, 'FaceColor', 'green')
xline(median(max_consec_channel_19(max_consec_channel_19 > 0)), 'k-', 'linewidth', 2)
xlabel('max consecutive hours channelized')
legend('treatment','median')
xlim([0 50])
ylim([0 0.15])
y_width=7.25;x_width=9.125;
set(gcf, 'PaperPosition', [0 0 x_width y_width]);
saveas(fig, '../figures/esurf_FigureB7.pdf')
