clear all; close all
%% The purpose of this script is to calculate channel depth and compensation timescale
% These data are plotted in Table 1 and Figure 3c in Sanks et al.
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

%% Calculate channel depth, aggradation, and compensation timescale 
% control 
C_R_18 = [];%empty array to be filled with Hc value along radial transects ever 5 mm from source 
C_Rv_18 = [];%empty array to be filled with std of Hc values along radial transects ever 5 mm from source 
Agg_18 = [];%empty array to be filled with aggradation rates (mm/hr) along radial transects ever 5 mm from source 
Aggv_18 = [];%empty array to be filled with variability of aggradation rates (mm/hr) along radial transects ever 5 mm from source 
TC_18 = [];%empty array to be filled with Tc value (hr) along radial transects ever 5 mm from source 
TCv_18 = [];%empty array to be filled with variability of Tc (hr) along radial transects ever 5 mm from source 

for i = 1:650;%loop to run through different radial distances from the end of the entrance channel.
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
    for j = 1:nanmax(size(xunit));
        xs_shot = ZD_18(xunit(j),yunit(j),:);
        is_shot = CM_18(xunit(j),yunit(j),:);
        zs18 = [zs18;xs_shot];
        is18 = [is18;is_shot];
    end
    zs18 = squeeze(zs18);%matrix of topo along radial transect
    is18 = squeeze(is18);%matrix of channel (yes/no) along radial transect
    zs2_18 = [];
    is2_18 = [];
    for j = 1:size(zs18,1);
        ztest = zs18(j,1);
        if ztest > 0
            zline = zs18(j,:);
            iline = is18(j,:);
            zs2_18 = [zs2_18;zline];
            is2_18 = [is2_18;iline];
        end
    end
    zs18 = zs2_18;
    is18 = is2_18;
    % Loop to find channel depths
    nw = 225;%window size
    Hc_list18 = [];
    Agg_list18 = [];
    Tc_list18 = [];
    for l = 1:1:size(zs18,2)-nw;
        zs_crop = zs18(:,l:l+nw-1);
        is_crop = is18(:,l:l+nw-1);
        Cdepth18 = [];
        for j = 1:size(zs_crop,2);
            iold = 0;
            for k = 2:size(zs_crop,1)-1;
                testi = is_crop(k,j);
                if testi > iold;%if, as going along transect, switch from non-channel to channel, start tracking topo
                    zc = zs_crop(k,j);
                end
                if testi == iold
                    if testi == 1;%if as going along transect, you stay channel node, keep tracking topo
                        zc = [zc;zs_crop(k,j)];
                    end
                end
                if testi < iold;%if, as going along transect, switch from channel to non-channel node, stop tracking topo and calculate Hc
                    cdepth = nanmax(zc) - nanmin(zc);
                    Cdepth18 = [Cdepth18;cdepth];
                end
                iold = testi;
            end
        end
        Hc_18 = prctile(Cdepth18,95);
        Hc_list18 = [Hc_list18;Hc_18];
        dz_18 = (nanmean(zs_crop(:,nw))-nanmean(zs_crop(:,1)))/nw;
        Agg_list18 = [Agg_list18;dz_18];
        Tc = Hc_18/dz_18;
        Tc_list18 = [Tc_list18;Tc];
    end
    Hc_18 = nanmean(Hc_list18);
    %if Hc_18 > 20;
    %    break
    %end
    Hc_v_18 = nanstd(Hc_list18);
    C_R_18 = [C_R_18;Hc_18];
    C_Rv_18 = [C_Rv_18;Hc_v_18];    
    dz_18 = nanmean(Agg_list18);
    dz_v_18 = nanstd(Agg_list18);
    Agg_18 = [Agg_18;dz_18];
    Aggv_18 = [Aggv_18;dz_v_18];
    tc_18 = nanmean(Tc_list18);
    tc_v_18 = nanstd(Tc_list18);
    TC_18 = [TC_18;tc_18];
    TCv_18 = [TCv_18;tc_v_18];
end

% treatment
C_R_19 = [];%empty array to be filled with Hc value along radial transects ever 5 mm from source 
C_Rv_19 = [];%empty array to be filled with std of Hc values along radial transects ever 5 mm from source 
Agg_19 = [];%empty array to be filled with aggradation rates (mm/hr) along radial transects ever 5 mm from source 
Aggv_19 = [];%empty array to be filled with variability of aggradation rates (mm/hr) along radial transects ever 5 mm from source 
TC_19 = [];%empty array to be filled with Tc value (hr) along radial transects ever 5 mm from source 
TCv_19 = [];%empty array to be filled with variability of Tc (hr) along radial transects ever 5 mm from source 

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
    for j = 1:nanmax(size(xunit))
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
    for j = 1:nanmax(size(xunit));
        xs_shot = ZD_19(xunit(j),yunit(j),:);
        is_shot = CM_19(xunit(j),yunit(j),:);
        zs19 = [zs19;xs_shot];
        is19 = [is19;is_shot];
    end
    zs19 = squeeze(zs19);%matrix of topo along radial transect
    is19 = squeeze(is19);%matrix of channel (yes/no) along radial transect
    zs2_19 = [];
    is2_19 = [];
    for j = 1:size(zs19,1);
        ztest = zs19(j,1);
        if ztest > 0
            zline = zs19(j,:);
            iline = is19(j,:);
            zs2_19 = [zs2_19;zline];
            is2_19 = [is2_19;iline];
        end
    end
    zs19 = zs2_19;
    is19 = is2_19;
    % Loop to find channel depths
    nw = 225;%window size
    Hc_list19 = [];
    Agg_list19 = [];
    Tc_list19 = [];
    for l = 1:1:size(zs19,2)-nw;
        zs_crop = zs19(:,l:l+nw-1);
        is_crop = is19(:,l:l+nw-1);
        Cdepth19 = [];
        for j = 1:size(zs_crop,2);
            iold = 0;
            for k = 2:size(zs_crop,1)-1;
                testi = is_crop(k,j);
                if testi > iold;%if, as going along transect, switch from non-channel to channel, start tracking topo
                    zc = zs_crop(k,j);
                end
                if testi == iold
                    if testi == 1;%if as going along transect, you stay channel node, keep tracking topo
                        zc = [zc;zs_crop(k,j)];
                    end
                end
                if testi < iold;%if, as going along transect, switch from channel to non-channel node, stop tracking topo and calculate Hc
                    cdepth = nanmax(zc) - nanmin(zc);
                    Cdepth19 = [Cdepth19;cdepth];
                end
                iold = testi;
            end
        end
        Hc_19 = prctile(Cdepth19,95);
        Hc_list19 = [Hc_list19;Hc_19];
        dz_19 = ((nanmean(zs_crop(:,nw))-nanmean(zs_crop(:,1)))/nw)/dt_19;
        Agg_list19 = [Agg_list19;dz_19];
        Tc = Hc_19/dz_19;
        Tc_list19 = [Tc_list19;Tc];
    end
    Hc_19 = nanmean(Hc_list19);
    Hc_v_19 = nanstd(Hc_list19);
    C_R_19 = [C_R_19;Hc_19];
    C_Rv_19 = [C_Rv_19;Hc_v_19];    
    dz_19 = nanmean(Agg_list19);
    dz_v_19 = nanstd(Agg_list19);
    Agg_19 = [Agg_19;dz_19];
    Aggv_19 = [Aggv_19;dz_v_19];
    tc_19 = nanmean(Tc_list19);
    tc_v_19 = nanstd(Tc_list19);
    TC_19 = [TC_19;tc_19];
    TCv_19 = [TCv_19;tc_v_19];
end

dist_18 = 1:650;
dist_19 = 1:650;
fig = figure();
subplot(2,1,1)
plot((dist_18*5)/1000, C_R_18, 'b-', 'linewidth',3);
grid on
grid minor
hold on
yline(nanmean(C_R_18), 'b-', 'linewidth', 1);
plot((dist_19*5)/1000, C_R_19, 'g-', 'linewidth',3);
yline(nanmean(C_R_19), 'g-', 'linewidth', 1);
xlim([0 3.5])
ylabel('channel depth (H_{c}; mm)');
legend('control', 'treatment', 'control mean', 'treatment mean');
subplot(2,1,2)
plot((dist_18*5)/1000, Agg_18, 'b-', 'linewidth',3);
grid on
grid minor
hold on
plot((dist_19*5)/1000, Agg_19, 'g-', 'linewidth',2);
yline((nanmean(Agg_18)), 'b-', 'linewidth',1);
yline((nanmean(Agg_19)), 'g-', 'linewidth',1);
xlabel('radial distance from the entrance channel (m)');
ylabel('channel aggradation rate (mm/hr)');
saveas(fig, 'chandepth_chanagg_dist.pdf');

figure()
plot(dist_18, TC_18, 'b-', 'linewidth',2);
grid on
grid minor
hold on
yline(nanmean(TC_18), 'b-', 'linewidth', 1);
plot(dist_19, TC_19, 'g-', 'linewidth',2);
yline(nanmean(TC_19), 'g-', 'linewidth', 1);
xlabel('radial distance from the entrance');
ylabel('compensation timescale (T_{c}; hours)');
legend('control', 'control mean', 'treatment', 'treatment mean');

%% Channel depth for Table 1
mean(C_R_18, 'omitnan')
std(C_R_18, 'omitnan')

mean(C_R_19, 'omitnan')
std(C_R_19, 'omitnan')