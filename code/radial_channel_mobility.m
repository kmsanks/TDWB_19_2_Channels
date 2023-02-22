clear all
%%Load data
cd '..\data'
load('CM_18.mat')
load('CM_19.mat')
load('ZD_18.mat')
load('ZD_19.mat')
cd '..\code'

%number of x locations on map
nx_18 = size(ZD_18, 1);
%number of y locations on map
ny_18 = size(ZD_18,2);
%number of time steps in data set
nt_18 = size(ZD_18,3);
dx = 5; %5 mm grid cells
dt_18 = 1; %delta t of time steps (hr)

%x grid node location of the entrance channel
xentrance_18 = 109;
%y grid node location of the entrance channel
yentrance_18 = 271;

%base level rise rate (mm/hr)
baselevel_rr = 0.25;
ocean_zero = 25; %ocean elevation at beginning of experiment (mm)

%replace empty channel maps in treatment with channel map during timestep
%after
for i = (560-1):-1:1; %I know 560 has a channel map, so I can start at 279 and replace with channel map that comes after, this will work for ones that have multiple no maps in a row
    if sum(sum(CM_19(:,:,i), 'omitnan')) == 0;
        CM_19(:,:,i) = CM_19(:,:,i+1);
    end
end

%number of x locations on map
nx_19 = size(CM_19,1);
%number of y locations on map
ny_19 = size(CM_19,2);
%number of time steps in data set
nt_19 = size(CM_19,3);
dt_19 = 1; %delta t of channel maps (hr)

%x grid node location of the entrance channel (x is down dip)
xentrance_19 = 214;
%y grid node location of the entrance channel (y is strike)
yentrance_19 = 397;

%%Create matrix of distances to channel entrance
[X18 Y18] = meshgrid(1:ny_18, 1:nx_18);
dd18 = sqrt((X18 - yentrance_18).^2 + (Y18 - xentrance_18).^2)*5;
%make everything outside of basin a NaN
tmp = zeros(796,522);
z = ZD_18(:,:,1);
z(z == 0.) = NaN;
tmp2 = z.*tmp;
tmp2(tmp2 == 0.) = 1;
dd18 = dd18.*tmp2;

[X19 Y19] = meshgrid(1:ny_19, 1:nx_19);
dd19 = sqrt((X19 - yentrance_19).^2 + (Y19 - xentrance_19).^2)*5;
%make everything outside of basin a NaN
tmp = zeros(750,747);
z = ZD_19(:,:,1);
tmp2 = z.*tmp;
tmp2(tmp2 == 0.) = 1;
dd19 = dd19.*tmp2;

%% Control Lateral Mobility
%look through 5mm radial distances and determine lateral mobility for each
%band

distance_18 = []; %radial distance from entrance channel
lat_mob_90_18 = []; %this will be mean hour it takes channel to visit 90% of radial band
mean_lat_mob_90_18 = [];
std_lat_mob_90_18 = [];
lat_mob_50_18 = []; %this will be mean hour it takes channel to visit 50% of radial band
mean_lat_mob_50_18 = [];
std_lat_mob_50_18 = [];
dist = 0:100:3000;
%%look to see if agg rate is divided by dt
for k = 1:length(0:100:3000)%loop to run through different radial distances from the end of the entrance channel.
    %section to find x,y nodes for radial transect and generate matrix of
    %cross section topo and channel (yes/no) data
    idx = dd18(dd18 >= dist(k) & dd18 < dist(k+1));
    radial_dd = dd18 >= dist(k) & dd18 < dist(k+1);
    area = length(idx)*2.5*10^-5; %area in m, not sure what this will be but fix here
    lat_chan_mob_18 = [];
    for i = 1:((size(CM_18,3)-60))
       for ii = i:size(CM_18,3)
           if ii == i
               chanMaps_n = CM_18(:,:,ii).*radial_dd;
               lat_chan_mob_18_i = 1-((sum(sum(chanMaps_n)*2.5*10^-5, 'omitnan'))/area);
               hour_18_i = (ii*2)-(i*2);
           else
               chanMaps_n = chanMaps_n + (CM_18(:,:,ii).*radial_dd);
               chanMaps_n(chanMaps_n > 1) = 1;
               lat_chan_mob_18_i = 1-(sum(sum(chanMaps_n)*2.5*10^-5, 'omitnan')/area);
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
    if length(mean_hr)>0
        hr = mean_hr(1);
        mean = mean_lat_chan_mob_18(hr,2);
        std = std_lat_chan_mob_18(hr,2);
        lat_mob_90_18(k) = hr;
        mean_lat_mob_90_18(k) = mean;
        std_lat_mob_90_18(k) = std;
    else 
        lat_mob_90_18(k) = NaN;
        mean_lat_mob_90_18(k) = NaN;
        std_lat_mob_90_18(k) = NaN;
    end
    %50% of area visited by channel
    [mean_hr] = find(mean_lat_chan_mob_18(:,2)<=0.5);
    if length(mean_hr)>0
        hr = mean_hr(1);
        mean = mean_lat_chan_mob_18(hr,2);
        std = std_lat_chan_mob_18(hr,2);
        lat_mob_50_18(k) = hr;
        mean_lat_mob_50_18(k) = mean;
        std_lat_mob_50_18(k) = std;
    else 
        lat_mob_50_18(k) = NaN;
        mean_lat_mob_50_18(k) = NaN;
        std_lat_mob_50_18(k) = NaN;
    end
    %distance from entrance channel
    distance_18(k) = dist(k)*0.001; %distance in m
end

%%Treatment Lateral Mobility
distance_19 = []; %radial distance from entrance channel
lat_mob_90_19 = []; %this will be mean hour it takes channel to visit 90% of radial band
mean_lat_mob_90_19 = [];
std_lat_mob_90_19 = [];
lat_mob_50_19 = []; %this will be mean hour it takes channel to visit 50% of radial band
mean_lat_mob_50_19 = [];
std_lat_mob_50_19 = [];
dist = 0:100:3500;
%%look to see if agg rate is divided by dt
for k = 1:length(0:100:3000)%loop to run through different radial distances from the end of the entrance channel.
    %section to find x,y nodes for radial transect and generate matrix of
    %cross section topo and channel (yes/no) data
    idx = dd19(dd19 >= dist(k) & dd19 < dist(k+1));
     %if k < 7
    radial_dd = dd19 >= dist(k) & dd19 < dist(k+1);
     %else
         %break
     %end
    area = length(idx)*2.5*10^-5; %area in m
    lat_chan_mob_19 = [];
    for i = 1:((size(CM_19,3)-60))
       for ii = i:size(CM_19,3)
           if ii == i
               chanMaps_n = CM_19(:,:,ii).*radial_dd;
               lat_chan_mob_19_i = 1-((sum(sum(chanMaps_n, 'omitnan')*2.5*10^-5))/area);
               hour_19_i = (ii*2)-(i*2);
           else
               chanMaps_n = chanMaps_n + (CM_19(:,:,ii).*radial_dd);
               chanMaps_n(chanMaps_n > 1) = 1;
               lat_chan_mob_19_i = 1-(sum(sum(chanMaps_n, 'omitnan')*2.5*10^-5)/area);
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
    if length(mean_hr)>0
        hr = mean_hr(1);
        mean = mean_lat_chan_mob_19(hr,2);
        std = std_lat_chan_mob_19(hr,2);
        lat_mob_90_19(k) = hr;
        mean_lat_mob_90_19(k) = mean;
        std_lat_mob_90_19(k) = std;
    else 
        lat_mob_90_19(k) = NaN;
        mean_lat_mob_90_19(k) = NaN;
        std_lat_mob_90_19(k) = NaN;
    end
    %50% of area visited by channel
    [mean_hr] = find(mean_lat_chan_mob_19(:,2)<=0.5);
    if length(mean_hr)>0
        hr = mean_hr(1);
        mean = mean_lat_chan_mob_19(hr,2);
        std = std_lat_chan_mob_19(hr,2);
        lat_mob_50_19(k) = hr;
        mean_lat_mob_50_19(k) = mean;
        std_lat_mob_50_19(k) = std;
    else 
        lat_mob_50_19(k) = NaN;
        mean_lat_mob_50_19(k) = NaN;
        std_lat_mob_50_19(k) = NaN;
    end
    %distance from entrance channel
    distance_19(k) = dist(k)*0.001; %distance in m
end

distance = dist*0.001;
distance = distance(2:end);
%%Lets plot the data now
figure
plot(distance_18, lat_mob_50_18, 'b-s')
hold on
plot(distance_18, lat_mob_90_18, 'c-^')
plot(distance_19, lat_mob_50_19, 'g-s')
plot(distance_19, lat_mob_90_19, 'g:^')
grid on
grid minor
legend('control 50%', 'control 90%', 'treatment 50%', 'treatment 90%')
xlabel('distance from the entrance channel (m)')
ylabel('lateral mobility (hrs)')

chanMaps_18_sum = [];
i = 1:50;
 for ii = 100:150
     if ii == 100
         chanMaps_n18 = CM_18(:,:,ii);
         %lat_chan_mob_19_i = 1-((nansum(nansum(chanMaps_n)*2.5*10^-5))/area);
         %hour_19_i = (ii*2)-(i*2);
     else
         chanMaps_n18 = chanMaps_n18 + CM_18(:,:,ii);
         chanMaps_n18(chanMaps_n18 > 1) = 1;
     end 
 end 