cd 'C:\Users\kmsanks\Box Sync\TDWB_19_2\TDB_18_1'
load('TDB_18_data.mat')
cd 'C:\Users\kmsanks\Box Sync\TDWB_19_2\Processed_Data\SurfaceProcesses\Topography\Timescales\LateralMobility'
%Channel maps and Z data
chanMaps_18 = C_maps;
chanMaps_18(:,:,298) = ~chanMaps_18(:,:,298);
chanMaps_18(:,:,333) = ~chanMaps_18(:,:,333);
ZD_18 = Z_maps(:,:,1:560);
clear('C_maps','B_maps','G_maps','R_maps','Z_maps');

%number of x locations on map
nx = size(ZD_18, 1);
%number of y locations on map
ny = size(ZD_18,2);
%number of time steps in data set
nt = size(ZD_18,3);
dx = 5; %5 mm grid cells
dt = 1; %delta t of time steps (hr)

%x grid node location of the entrance channel
xentrance = 109;
%y grid node location of the entrance channel
yentrance = 271;

%base level rise rate (mm/hr)
baselevel_rr = 0.25;
ocean_zero = 25; %ocean elevation at beginning of experiment (mm)

%fill empty matrics
pland_18 = [];
t = [];
for i = 1:nt;
    t = [t;i*dt];
    z = ZD_18(:,:,i);
    z = z - (baselevel_rr*i*dt+ocean_zero);%
    z(z > 0) = 1;
    z(z < 0) = 0;
    pland_18(:,:,i) = z;
end 

load('Z:\TDWB_19_2\Processed_Data\TDWB_19_2_lidar_scans\ZD_19_2_dry.mat'); %load topography array. In here should be a 3D topo array called ZD, oriented space x space x time
load('C:\Users\kmsanks\Box Sync\TDWB_19_2\Processed_Data\SurfaceProcesses\Images\ChannelMaps\TDWB_19_2_chanMaps.mat');
cd 'C:\Users\kmsanks\Box Sync\TDWB_19_2\Processed_Data\SurfaceProcesses\Topography\Timescales\LateralMobility'

chanMaps_19 = TDWB_19_2_chanMaps;
clear('TDWB_19_2_chanMaps');
%number of x locations on map
nx = size(ZD_19_2_dry,1);
%number of y locations on map
ny = size(ZD_19_2_dry,2);
%number of time steps in data set
nt = size(ZD_19_2_dry,3);
dx = 5; %5 mm grid cells
dt = 2; %delta t of time steps (hr)
%Set the boundary conditions
baselevel_rr = 0.25; %base level rise rate (mm/hr)
ocean_zero = 25; %ocean elevation at beginning of experiment (mm)
%loop through all x locations on the delta top
pland_19 = [];
t_19 = [];
for i = 1:nt;
    t_19 = [t;i*dt]; %time
    z = ZD_19_2_dry(:,:,i); %elevations
    z = z - (baselevel_rr*i*dt+ocean_zero); %ele
    z(z > 0) = 1;
    z(z < 0) = 0;
    pland_19(:,:,i) = z;
end 
%add all cells through time
pland_18_sum = sum(pland_18,3);
pland_18_sum(pland_18_sum == 0) = NaN;
pland_18_frac = pland_18_sum/565;
frac_18 = [];
area_pland_18 = [];
for i = round(linspace(0,1,101),2);
    frac_18 = [frac_18;i];
    area = pland_18_frac(pland_18_frac >= i); %& pland_18_frac < i+1
    area(area > 0) = 1;
    areasum = sum(area);
    area_pland_18_i = (areasum)*2.5*10^-5; %area of 1 pixel
    area_pland_18 = [area_pland_18; area_pland_18_i];     
end 

pland_19_sum = sum(pland_19,3);
pland_19_sum(pland_19_sum == 0) = NaN;
pland_19_frac = pland_19_sum/281;
frac_19 = [];
area_pland_19 = [];
for i = round(linspace(0,1,101),2);
    frac_19 = [frac_19;i];
    area = pland_19_frac(pland_19_frac >= i); %& pland_18_frac < i+1
    area(area > 0) = 1;
    areasum = sum(area);
    area_pland_19_i = (areasum)*2.5*10^-5; %area of 1 pixel
    area_pland_19 = [area_pland_19; area_pland_19_i];     
end 

frac_area_array = [frac_18, area_pland_18];
%50 percent area
frac50_index = find(frac_area_array(:,1) == [0.50]);
area_frac50_18 = frac_area_array(frac50_index,2);

frac_area_array_19 = [frac_19, area_pland_19];
%50 percent area
frac50_index_19 = find(frac_area_array_19(:,1) == [0.50]);
area_frac50_19 = frac_area_array_19(frac50_index_19,2);

%%Now lets look at channel mobility
%We want to know (starting at each hour) how long it takes to visit 95% of
%the delta top
pland_18_frac_nan = pland_18_frac;
pland_18_frac_nan(pland_18_frac_nan < 0.5) = NaN;
pland_18_frac_nan(pland_18_frac_nan > 0) = 1;
chanMaps_18_crop = [];
for i = 1:((size(chanMaps_18,3)));
    chanMaps_18_nan = chanMaps_18(:,:,i).*pland_18_frac_nan;
    chanMaps_18_crop(:,:,i) = chanMaps_18_nan;
end
% %remove every other element to make the same timestep as 19
% chanMaps_18_lat = chanMaps_18_crop(:,:,1:2:end);

lat_chan_mob_18 = [];
time_2_66_18 = [];
for i = 1:((size(chanMaps_18_crop,3)-60));
    tmp_mob = [];
    for ii = i:size(chanMaps_18_crop,3);
        area = area_frac50_18; %0.5 land above sea level 2.1811 m^2
        if ii == i;
            chanMaps_n = chanMaps_18_crop(:,:,ii);
            lat_chan_mob_18_i = 1-((nansum(nansum(chanMaps_n)*2.5*10^-5))/area);
            hour_18_i = (ii*2)-(i*2);
        else
            chanMaps_n = chanMaps_n + chanMaps_18_crop(:,:,ii);
            chanMaps_n(chanMaps_n > 1) = 1;
            lat_chan_mob_18_i = 1-(nansum(nansum(chanMaps_n)*2.5*10^-5)/area);
            hour_18_i = ii-i;
            lat_chan_mob_18_i = [hour_18_i, lat_chan_mob_18_i];
            lat_chan_mob_18 = [lat_chan_mob_18;lat_chan_mob_18_i];
            tmp_mob = [tmp_mob;lat_chan_mob_18_i];
        end
    end
    [idx] = tmp_mob(tmp_mob(:,2)<= 0.34);
    if isempty(idx)
        time_2_66_18(i) = NaN;
    else 
        time_2_66_18(i) = min(idx);
    end
    %time_2_95 = lat_chan_mob_18(idx,1);
end 

%fraction unmodified
%TDB-19
%%Replace timesteps with no channel maps with the channel map from the
%%next time step for the treatment experiment
for i = (560-1):-1:1; %I know 560 has a channel map, so I can start at 279 and replace with channel map that comes after, this will work for ones that have multiple no maps in a row
    if nansum(nansum(chanMaps_19(:,:,i))) == 0;
        chanMaps_19(:,:,i) = chanMaps_19(:,:,i+1);
    end
end
pland_19_frac_nan = pland_19_frac;
pland_19_frac_nan(pland_19_frac_nan < 0.5) = NaN;
pland_19_frac_nan(pland_19_frac_nan > 0) = 1;
chanMaps_19_crop = [];
for i = 1:((size(chanMaps_19,3)));
    chanMaps_19_nan = chanMaps_19(:,:,i).*pland_19_frac_nan;
    chanMaps_19_crop(:,:,i) = chanMaps_19_nan;
end

lat_chan_mob_19 = [];
time_2_67_19 = [];
for i = 1:((size(chanMaps_19_crop,3)-60));
    tmp_mob = [];
    for ii = i:size(chanMaps_19_crop,3);
        area = area_frac50_19; %0.5 land above sea level 1.9329 m^2
        if ii == i;
            chanMaps_n = chanMaps_19_crop(:,:,ii);
            lat_chan_mob_19_i = 1-((nansum(nansum(chanMaps_n)*2.5*10^-5))/area);
            hour_19_i = ii-i;
        else
            chanMaps_n = chanMaps_n + chanMaps_19_crop(:,:,ii);
            chanMaps_n(chanMaps_n > 1) = 1;
            lat_chan_mob_19_i = 1-(nansum(nansum(chanMaps_n)*2.5*10^-5)/area);
            hour_19_i = ii-i;
            lat_chan_mob_19_i = [hour_19_i, lat_chan_mob_19_i];
            lat_chan_mob_19 = [lat_chan_mob_19;lat_chan_mob_19_i];
            tmp_mob = [tmp_mob;lat_chan_mob_19_i];
        end
    end
    [idx] = tmp_mob(tmp_mob(:,2)<= 0.33);
    if isempty(idx)
        time_2_67_19(i) = NaN;
    else 
        time_2_67_19(i) = min(idx);
    end 
end

fig = figure();
plot(1:length(time_2_66_18), time_2_66_18, 'bo')
hold on
plot(1:length(time_2_67_19), time_2_67_19, 'go')
xlabel('begin hour')
ylabel('time for channel to visit 1 e-fold area of the delta top (hrs)')
legend('control (66%)', 'treatment (67%)')
set(gcf, 'PaperUnits', 'inches');
y_width=7.25;x_width=9.125
set(gcf, 'PaperPosition', [0 0 x_width y_width]);
saveas(fig, 'time_2_66_latmob.pdf')    
    
[ud,ix,iy] = unique(lat_chan_mob_18(:,1));  
mean_lat_chan_mob_18 = [ud, accumarray(iy,lat_chan_mob_18(:,2),[],@mean)];
std_lat_chan_mob_18 = [ud, accumarray(iy,lat_chan_mob_18(:,2),[],@std)];

[ud,ix,iy] = unique(lat_chan_mob_19(:,1));  
mean_lat_chan_mob_19 = [ud, accumarray(iy,lat_chan_mob_19(:,2),[],@mean)];
std_lat_chan_mob_19 = [ud, accumarray(iy,lat_chan_mob_19(:,2),[],@std)];
std_lat_chan_mob_19(std_lat_chan_mob_19 < 0) = 0;

figure()
plot(mean_lat_chan_mob_18(:,1), mean_lat_chan_mob_18(:,2), 'bo');
hold on
plot(mean_lat_chan_mob_19(:,1), mean_lat_chan_mob_19(:,2), 'go');
grid on
grid minor
xlim([0 350]);
xlabel('measurement window (hr)');
ylabel('fraction visited  by channel (-)');
legend('control', 'treatment');

%Log-space
figure()
semilogy(mean_lat_chan_mob_18(:,1), mean_lat_chan_mob_18(:,2), 'bo');
hold on
semilogy(mean_lat_chan_mob_19(:,1), mean_lat_chan_mob_19(:,2), 'go');
grid on
grid minor
xlim([0 560]);
xlabel('measurement window (hr)');
ylabel('log(fraction visited by channel (-))');
legend('control', 'treatment');

%x grid node location of the entrance channel
x18entrance = 109;
%y grid node location of the entrance channel
y18entrance = 271;
%x grid node location of the entrance channel (x is down dip)
x19entrance = 214;
%y grid node location of the entrance channel (y is strike)
y19entrance = 397;

%differnce to plot contours on each other
xdif = x19entrance - x18entrance;
ydif = y19entrance - y18entrance;

%subtract xdif and ydif from 19 matrix
pland_19_frac_shift = pland_19_frac(xdif:end, ydif:end, :);
pland_19_frac_nan_shift = pland_19_frac_nan(xdif:end, ydif:end, :);
%plot contours of 18 and 18 50% land
v = [0.5,0.5];
figure()
imagesc(pland_18_frac_nan)
contour(pland_18_frac,v, 'b-')
hold on
%imagesc(pland_19_frac_nan_shift)
%alpha(pland_19_frac_nan_shift, 'clear', 0.1)
contour(pland_19_frac_shift, v, 'g-')
legend('TDB-18-1', 'TDWB-19-2')
set(gca, 'Ydir', 'reverse');

%fraction unmodified
%crop elevations to area mask
ZD_18_crop = ZD_18(:,:,1:2:end).*pland_18_frac_nan;
ZD_19_crop = ZD_19_2_dry.*pland_19_frac_nan;

nt=280;
f_um_18 = [];
for i = 1:nt-1;
    for ii = i:nt;
        area = area_frac50_18; %0.5 land above sea level 1.9329 m^2
        dz_18_i = ZD_18_crop(:,:,ii)-ZD_18_crop(:,:,i);
        dz_18_i(dz_18_i < 1) = 0;
        dz_18_i(dz_18_i > 1) = 1;
        f_um_18_i = 1-(nansum(nansum(dz_18_i)*2.5*10^-5)/area);
        ts_i = ii-i;
        f_um_18_i = [ts_i, f_um_18_i];
        f_um_18 = [f_um_18;f_um_18_i];
    end
end

nt=281;
f_um_19 = [];
for i = 1:(size(ZD_19_crop,3)-1);
    for ii = i:size(ZD_19_crop,3);
        area = area_frac50_19; %0.5 land above sea level 1.9329 m^2
        dz_19_i = ZD_19_crop(:,:,ii)-ZD_19_crop(:,:,i);
        dz_19_i(dz_19_i < 1) = 0;
        dz_19_i(dz_19_i > 1) = 1;
        f_um_19_i = 1-(nansum(nansum(dz_19_i)*2.5*10^-5)/area);
        ts_i = ii-i;
        f_um_19_i = [ts_i, f_um_19_i];
        f_um_19 = [f_um_19;f_um_19_i];
    end
end

[ud,ix,iy] = unique(f_um_18(:,1));  
mean_f_um_18 = [ud, accumarray(iy,f_um_18(:,2),[],@mean)];
std_f_um_18 = [ud, accumarray(iy,f_um_18(:,2),[],@std)];

[ud,ix,iy] = unique(f_um_19(:,1));  
mean_f_um_19 = [ud, accumarray(iy,f_um_19(:,2),[],@mean)];
std_f_um_19 = [ud, accumarray(iy,f_um_19(:,2),[],@std)];

figure()
plot(mean_f_um_18(:,1), mean_f_um_18(:,2), 'bo');
hold on
plot(mean_f_um_19(:,1), mean_f_um_19(:,2), 'go');
grid on
grid minor
xlim([0 150]);
xlabel('measurement window (hr)');
ylabel('fraction unmodified (-)');
legend('control', 'treatment');


%plot with standard deviation
%lateral channel mobility
array18_lat = [(1:559); mean_lat_chan_mob_18(:,2)'; std_lat_chan_mob_18(:,2)'];
cols = any(isnan(array18_lat),1);
array18_lat(:,cols) = [];

array19_lat = [(1:559);mean_lat_chan_mob_19(:,2)'; std_lat_chan_mob_19(:,2)'];
cols = any(isnan(array19_lat),1);
array19_lat(:,cols) = [];

%fill standard deviation
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

%fraction unmodified
array18_fum = [(1:280); mean_f_um_18(:,2)'; std_f_um_18(:,2)'];
cols = any(isnan(array18_fum),1);
array18_fum(:,cols) = [];

array19_fum = [(1:281);mean_f_um_19(:,2)'; std_f_um_19(:,2)'];
cols = any(isnan(array19_fum),1);
array19_fum(:,cols) = [];

%fill standard deviation
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

fig = figure()
%plot lateral channel mobility
subplot(2,2,1)
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
grid on
grid minor
ylabel('fraction not visited by channel (-)')
xlabel('measurement window (hr)') 
legend('control mean', 'treatment mean', 'control stdev', 'treatment stdev')
%plot fraction unmodified 
subplot(2,2,2)
plot(x18_fum, y18_fum, 'b', 'LineWidth', 2)
hold on
plot(x19_fum, y19_fum, 'g', 'LineWidth', 2)
patch([x18_fum fliplr(x18_fum)], [curve1_18_fum fliplr(curve2_18_fum)], 'b')
patch([x19_fum fliplr(x19_fum)], [curve1_19_fum fliplr(curve2_19_fum)], 'g')
alpha(0.15)
plot(x18_fum, y18_fum, 'b', 'LineWidth', 2)
plot(x19_fum, y19_fum, 'g', 'LineWidth', 2)
ylim([0 1])
xlim([0 100])
grid on
grid minor
ylabel('fraction unmodified (-)')
xlabel('measurment window (hr)') 
legend('control mean', 'treatment mean', 'control stdev', 'treatment stdev')
%plot in semilogy space
subplot(2,2,3)
semilogy(mean_lat_chan_mob_18(:,1), mean_lat_chan_mob_18(:,2), 'bo');
hold on
semilogy(mean_lat_chan_mob_19(:,1), mean_lat_chan_mob_19(:,2), 'go');
grid on
grid minor
xlim([0 560]);
xlabel('measurement window (hr)');
ylabel('log(fraction not visited by channel (-))');
legend('control', 'treatment');
subplot(2,2,4)
semilogy(mean_f_um_18(:,1), mean_f_um_18(:,2), 'bo');
grid on
grid minor
hold on
semilogy(mean_f_um_19(:,1), mean_f_um_19(:,2), 'go');
grid on
grid minor
xlim([0 150]);
xlabel('measurement window (hr)');
ylabel('log(fraction unmodified (-))');
legend('control', 'treatment');




%patch([min(xlim) max(xlim) max(xlim) min(xlim)], [ybars(1) ybars(1), ybars(2) ybars(2)], [0.8 0.8 0.8])
set(gcf, 'PaperUnits', 'inches');
y_width=7.25;x_width=9.125
set(gcf, 'PaperPosition', [0 0 x_width y_width]);
saveas(fig, 'fraction_unmodified.pdf')




%exponential fit for e-folding ts
x = mean_f_um_18(:,1);
x1 = x(1:60);
y = mean_f_um_18(:,2);
y1 = y(1:60);
f0 = ezfit(x1, y1, 'exp');
xx = linspace(2,61,30);
figure()
semilogy(x1,y1,'o',xx,f0(xx),'r-');

efold18 = -1/f0.m(2);

%exponential fit for e-folding ts
x = mean_f_um_19(:,1);
x2 = x(1:80);
y = mean_f_um_19(:,2);
y2 = y(1:80);
f1 = ezfit(x2, y2, 'exp');
xx1 = round(linspace(1,80,40));
figure()
semilogy(x2,y2,'o',xx1,f1(xx1),'g-');

efold19 = -1/f1.m(2);

%exponential fit for lateral channel mobility
x = mean_lat_chan_mob_18(:,1);
x1 = x(1:380);
y = mean_lat_chan_mob_18(:,2);
y1 = y(1:380);
f0 = ezfit(x1, y1, 'exp');
xx = linspace(2,61,30);
figure()
semilogy(x1,y1,'o',xx,f0(xx),'r-');

latfold18 = -1/f0.m(2);

%exponential fit for lateral mobility
x = mean_lat_chan_mob_19(:,1);
x2 = x(1:240);
y = mean_lat_chan_mob_19(:,2);
y2 = y(1:240);
f1 = ezfit(x2, y2, 'exp');
xx1 = round(linspace(1,80,40));
figure()
semilogy(x2,y2,'o',xx1,f1(xx1),'g-');

latfold19 = -1/f1.m(2);
%Log-space
figure()
semilogy(mean_f_um_18(:,1), mean_f_um_18(:,2), 'bo');
hold on
semilogy(mean_f_um_19(:,1), mean_f_um_19(:,2), 'go');
grid on
grid minor
xlim([0 150]);
xlabel('measurement window (hr)');
ylabel('log(fraction unmodified (-))');
legend('control', 'treatment');








%exponential fit for e-folding ts
x = mean_lat_chan_mob_18(:,1);
y = mean_lat_chan_mob_18(:,2);
g = fittype('a-b*exp(-c*x)');
f0 = fit(x, y, g,'StartPoint',[[ones(size(x)), -exp(-x)]\y; 1]);
xx = linspace(1,559,1);
figure()
plot(x,y,'o',xx,f0(xx),'r-');

%exponential fit for e-folding ts
x1 = mean_lat_chan_mob_19(:,1);
y1 = mean_lat_chan_mob_19(:,2);
f1 = fit(x1, y1, g,'StartPoint',[[ones(size(x1)), -exp(-x1)]\y1; 1]);
xx1 = linspace(1,559,1);
figure()
plot(x1,y1,'o',xx1,f1(xx1),'g-');

%channel density and histogram
chan_dens18 = nansum(chanMaps_18_crop, 3);
chan_dens19 = nansum(chanMaps_19_crop, 3);

fig = figure
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
y_width=7.25;x_width=9.125
set(gcf, 'PaperPosition', [0 0 x_width y_width]);
saveas(fig, 'channel_density.pdf')


