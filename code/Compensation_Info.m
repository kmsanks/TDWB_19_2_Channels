clear all; close all
%% Definitions
blrise = 0.25;%base level rise rate (mm/hr)
dx = 5;%cell spacing in x direction (mm)
dy = 5;%cell spacing in y direction (mm)
dt = 2;%time between scans (hr)

%% Control Def
cd '../data'
load('ZD_18.mat')
ZD_18 = ZD_18(:,:,1:2:560);%3D array of topography (x,y,t array with values in z)
cd '../code'
z = ZD_18;%renaming ZD array
clear ZD_18
nx = size(z,1);%number of x nodes
ny = size(z,2);%number of y nodes
nt = size(z,3);%number of t nodes
xcenter = 109;%x node of end of entrance channel
ycenter = 271;%y node of end of entrance channel

%% Main loop control
Tc18 = [];%empty array to be filled with Tc value along radial transects ever 5 mm from source 
KAPPA_ST18 = [];%empty array to be filled with short-term compensation index ever 5 mm from source 
KAPPA_LT18 = [];%empty array to be filled with long-term compensation index ever 5 mm from source 
Agg18 = [];%empty array to be filled with aggradation rates (mm/hr) along radial transects ever 5 mm from source 
for i = 1:337%loop to run through different radial distances from the end of the entrance channel.
    i
    %%section to find x,y nodes for radial transect and generate matrix of
    %%cross section topo and channel (yes/no) data
    th = 0:1/i:2*pi;
    xunit = i * cos(th) + xcenter;
    yunit = i * sin(th) + ycenter;
    xunit = round(xunit);
    yunit = round(yunit);
    xunit2 = [];
    yunit2 = [];
    for j = 1:max(size(xunit))
        xc = xunit(j);
        yc = yunit(j);
        if xc >= 1;
            if yc >= 1;
                if xc <= nx
                    if yc <= ny
                        xunit2 = [xunit2;xc];
                        yunit2 = [yunit2;yc];
                    end
                end
            end
        end
    end
    xunit = xunit2;
    yunit = yunit2;
    XY = [yunit xunit];
    XY = sortrows(XY);
    xunit = XY(:,2);
    yunit = XY(:,1);
    zs = [];
    for j = 1:max(size(xunit));
        xs_shot = z(xunit(j),yunit(j),:);
        zs = [zs;xs_shot];
    end
    zs = squeeze(zs);%matrix of topo along radial transect
    zs2 = [];
    for j = 1:size(zs,1);
        ztest = zs(j,1);
        if ztest > 0
            zline = zs(j,:);
            zs2 = [zs2;zline];
        end
    end
    zs = zs2;%space x time
    %% Loop to find aggradation rate
    nw = nt-1;%window size, was 225
    Agg_list = [];
    for l = 1:1:size(zs,2)-nw;
        zs_crop = zs(:,l:l+nw-1);
        dz = (mean(zs_crop(:,nw))-mean(zs_crop(:,1)))/nw;
        Agg_list = [Agg_list;dz];
    end
    dz = mean(Agg_list);
    Agg18 = [Agg18;dz];
    [t,sss,TC,kappa_st,kappa_lt] = sigma_ss(zs,dt,blrise);
    Tc18 = [Tc18;TC];
    KAPPA_ST18 = [KAPPA_ST18;kappa_st];
    KAPPA_LT18 = [KAPPA_LT18;kappa_lt];
    hold off
end

%% Treatment def
cd '../data'
load('ZD_19.mat')
cd '../code'
z = ZD_19;%renaming ZD array
clear ZD_19
nx = size(z,1);%number of x nodes
ny = size(z,2);%number of y nodes
nt = size(z,3);%number of t nodes
xcenter = 214;%x node of end of entrance channel
ycenter = 397;%y node of end of entrance channel

%% Main loop
Tc19 = [];%empty array to be filled with Tc value along radial transects ever 5 mm from source 
KAPPA_ST19 = [];%empty array to be filled with short-term compensation index ever 5 mm from source 
KAPPA_LT19 = [];%empty array to be filled with long-term compensation index ever 5 mm from source 
Agg19 = [];%empty array to be filled with aggradation rates (mm/hr) along radial transects ever 5 mm from source 
for i = 1:337%loop to run through different radial distances from the end of the entrance channel.
    i
    %%section to find x,y nodes for radial transect and generate matrix of
    %%cross section topo and channel (yes/no) data
    th = 0:1/i:2*pi;
    xunit = i * cos(th) + xcenter;
    yunit = i * sin(th) + ycenter;
    xunit = round(xunit);
    yunit = round(yunit);
    xunit2 = [];
    yunit2 = [];
    for j = 1:max(size(xunit))
        xc = xunit(j);
        yc = yunit(j);
        if xc >= 1;
            if yc >= 1;
                if xc <= nx
                    if yc <= ny
                        xunit2 = [xunit2;xc];
                        yunit2 = [yunit2;yc];
                    end
                end
            end
        end
    end
    xunit = xunit2;
    yunit = yunit2;
    XY = [yunit xunit];
    XY = sortrows(XY);
    xunit = XY(:,2);
    yunit = XY(:,1);
    zs = [];
    for j = 1:max(size(xunit));
        xs_shot = z(xunit(j),yunit(j),:);
        zs = [zs;xs_shot];
    end
    zs = squeeze(zs);%matrix of topo along radial transect
    zs2 = [];
    for j = 1:size(zs,1);
        ztest = zs(j,1);
        if ztest > 0
            zline = zs(j,:);
            zs2 = [zs2;zline];
        end
    end
    zs = zs2;%space x time
    %% Loop to find aggradation rate
    nw = nt-1;%window size, was 225
    Agg_list = [];
    for l = 1:1:size(zs,2)-nw;
        zs_crop = zs(:,l:l+nw-1);
        dz = (mean(zs_crop(:,nw))-mean(zs_crop(:,1)))/nw;
        Agg_list = [Agg_list;dz];
    end
    dz = mean(Agg_list);
    Agg19 = [Agg19;dz];
    [t,sss,TC,kappa_st,kappa_lt] = sigma_ss(zs,dt,blrise);
    Tc19 = [Tc19;TC];
    KAPPA_ST19 = [KAPPA_ST19;kappa_st];
    KAPPA_LT19 = [KAPPA_LT19;kappa_lt];
    hold off
end


%% Figures
dist = ((1:1:337)*.005)';
dist18 = ((1:1:336)*.005)';
figure()
plot(dist18, Tc18, 'b-', 'LineWidth', 2)
hold on 
plot(dist, Tc19, 'g-', 'LineWidth', 2)
ylabel('compensation timescale (hrs)')
xlabel('distance from apex (m)')
legend('control', 'treatment')

figure()
plot(dist, Agg18, 'b-', 'LineWidth', 2)
hold on
plot(dist, Agg19, 'g-', 'LineWidth', 2)
ylabel('aggradation rate (mm/2-hr)')
xlabel('distance from apex (m)')
legend('control', 'treatment')

figure()
plot(dist18, KAPPA_ST18, 'b:', 'LineWidth', 2)
hold on
plot(dist, KAPPA_ST19, 'g:', 'LineWidth', 2)
plot(dist18, KAPPA_LT18, 'b-', 'LineWidth', 2)
plot(dist, KAPPA_LT19, 'g-', 'LineWidth', 2)
xlabel('distance from apex (m)')
ylabel('kappa (-)')
legend('control short-term', 'treatment short-term',...
    'control long-term', 'treatment long-term')