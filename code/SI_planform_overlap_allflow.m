clear all; close all;

%% Load data for the whole program.
cd 'C:\Users\kmsanks\Box Sync\TDWB_19_2\TDB_18_1\'
load('TDB_18_data.mat');
load('ZD.mat');
cd 'C:\Users\kmsanks\Box Sync\TDWB_19_2\Processed_Data\SurfaceProcesses\Topography\Timescales\LateralMobility'

ZD_18 = ZD(:,:,1:560);
clear('C_maps','B_maps','G_maps','R_maps','Z_maps', 'ZD');

%Load treatment data
load('Z:\TDWB_19_2\Processed_Data\TDWB_19_2_lidar_scans\ZD_19_2_dry.mat'); %load topography array. In here should be a 3D topo array called ZD, oriented space x space x time
cd 'C:\Users\kmsanks\Box Sync\TDWB_19_2\Processed_Data\SurfaceProcesses\Topography\Timescales\LateralMobility'

ZD_19 = ZD_19_2_dry;
clear('ZD_19_2_dry');

%load flow data
load('C:\Users\kmsanks\Box Sync\TDWB_19_2\Processed_Data\TDWB_19_2_Matrices\flowscreen18.mat');
load('C:\Users\kmsanks\Box Sync\TDWB_19_2\Processed_Data\TDWB_19_2_Matrices\flowscreen19.mat');
flowscreen18 = flowscreen18(:,:,1:560);

%%Calulate active floodplain area and number of wet and pixels for each time step
dt_18 = 1;
dt_19 = 2;

nt_18 = size(ZD_18,3);
nt_19 = size(ZD_19,3);

baselevel_rr = 0.25;
ocean_zero = 25;

pland_18 = [];
t_18 = [];
for i = 1:nt_18;
    t_18 = [t_18;i*dt_18];
    z = ZD_18(:,:,i);
    z = z - (baselevel_rr*i*dt_18+ocean_zero);%
    z(z > 0) = 1;
    z(z < 0) = 0;
    pland_18(:,:,i) = z;
end 

pland_19 = [];
t_19 = [];
for i = 1:nt_19;
    t_19 = [t_19;i*dt_19]; %time
    z = ZD_19(:,:,i); %elevations
    z = z - (baselevel_rr*(i-1)*dt_19+ocean_zero); %ele
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

pland_18_frac_nan = pland_18_frac;
pland_18_frac_nan(pland_18_frac_nan < 0.5) = NaN;
pland_18_frac_nan(pland_18_frac_nan > 0) = 1;
flowscreen_18_crop = [];
for i = 1:((size(flowscreen18,3)));
    flowscreen_18_nan = flowscreen18(:,:,i).*pland_18_frac_nan;
    flowscreen_18_crop(:,:,i) = flowscreen_18_nan;
end

pland_19_frac_nan = pland_19_frac;
pland_19_frac_nan(pland_19_frac_nan < 0.5) = NaN;
pland_19_frac_nan(pland_19_frac_nan > 0) = 1;
flowscreen_19_crop = [];
for i = 1:((size(flowscreen19,3)));
    flowscreen_19_nan = flowscreen19(:,:,i).*pland_19_frac_nan;
    flowscreen_19_crop(:,:,i) = flowscreen_19_nan;
end

%we need to make flow 1 and non-flow 0
flowscreen_18_mat = 1 - flowscreen_18_crop;
flowscreen_19_mat = 1 - flowscreen_19_crop;

%Total area
A_18 = nansum(nansum(pland_18_frac_nan));
A_19 = nansum(nansum(pland_19_frac_nan));

%some of the 19 flow maps do not have flow
for i = 1:size(flowscreen_19_mat,3)
    if nansum(nansum(flowscreen_19_mat(:,:,i))) == 0
        flowscreen_19_mat(:,:,i) = NaN;
    end
end

%Wet and dry pixels
wetA_18 = [];
dryA_18 = [];
for i = 1:size(flowscreen_18_mat,3)
    wetA = nansum(nansum(flowscreen_18_mat(:,:,i)));
    dryA = A_18-wetA;
    wetA_18 = [wetA_18;wetA];
    dryA_18 = [dryA_18;dryA];
end

%Wet and dry pixels
wetA_19 = [];
dryA_19 = [];
for i = 1:size(flowscreen_19_mat,3)
    wetA = nansum(nansum(flowscreen_19_mat(:,:,i)));
    dryA = A_19-wetA;
    wetA_19 = [wetA_19;wetA];
    dryA_19 = [dryA_19;dryA];
end

%%Now lets calculate the planform overlap
% Load each baseline.
changedpix18 = [];
hour18 = [];
Phi18 = [];
corr18 = [];
for z=1:size(flowscreen_18_mat,3)
    baseline = flowscreen_18_mat(:,:,z);   
    for f=z:size(flowscreen_18_mat,3)
% STEP 1: Find the number of decorrelated pixels.
        transient = flowscreen_18_mat(:,:,f);
        changedpix = nansum(nansum(abs(baseline-transient)));
        hour = f-z;
        changedpix18 = [changedpix18;changedpix];
        hour18 = [hour18;hour];
        % Note: In the channel-maps, 0=active floodplain and 1=channel. So
        % in changedpix, 0=no change and 1=change. 
        % In order to analyze the data, a single column matrix is made,
        % containing each of the sums.
        %npix_decor_outmat(f,2)=changedpix;
% STEP 2: Build scaling parameters and apply the overlap.
        Phi = (wetA_18(z)/A_18*dryA_18(f)/A_18)+(dryA_18(z)/A_18)*wetA_18(f)/A_18; 
        corr = 1-changedpix/(Phi*A_18); %A is just area of active floodplain (above rsl?)
        Phi18 = [Phi18;Phi];
        corr18 = [corr18;corr];
        % STEP 3: Set the baseline time-step to occur at t=0.
%         corr_zero_outmat(f,2)=corr_outmat(f,2);
        
    end
end

changedpix19 = [];
hour19 = [];
Phi19 = [];
corr19 = [];
for z=1:size(flowscreen_19_mat,3)
    baseline = flowscreen_19_mat(:,:,z);   
    for f=z:size(flowscreen_19_mat,3)
% STEP 1: Find the number of decorrelated pixels.
        transient = flowscreen_19_mat(:,:,f);
        if nansum(nansum(transient)) == 0 | nansum(nansum(baseline)) == 0
           changedpix = NaN;
           Phi = NaN;
           corr = NaN;
        else
        changedpix = nansum(nansum(abs(baseline-transient)));
        Phi = (wetA_19(z)/A_19*dryA_19(f)/A_19)+(dryA_19(z)/A_19)*wetA_19(f)/A_19; 
        corr = 1-changedpix/(Phi*A_19); %A is just area of active floodplain (above rsl?)
        end 
        hour = f-z;
        changedpix19 = [changedpix19;changedpix];
        hour19 = [hour19;hour];
        % Note: In the channel-maps, 0=active floodplain and 1=channel. So
        % in changedpix, 0=no change and 1=change. 
        % In order to analyze the data, a single column matrix is made,
        % containing each of the sums.
        %npix_decor_outmat(f,2)=changedpix;
% STEP 2: Build scaling parameters and apply the overlap.
        
        Phi19 = [Phi19;Phi];
        corr19 = [corr19;corr];
        % STEP 3: Set the baseline time-step to occur at t=0.
%         corr_zero_outmat(f,2)=corr_outmat(f,2);
        
    end
end

%changedpix19(changedpix19 == 0) = NaN;
%Phi19(Phi

planform_overlap_ob_18 = [hour18,changedpix18,Phi18,corr18];
planform_overlap_ob_19 = [hour19,changedpix19,Phi19,corr19];

%Control Data
[ud,ix,iy] = unique(planform_overlap_ob_18(:,1));  
%changed pix
mean_changedpix18 = [ud, accumarray(iy,planform_overlap_ob_18(:,2),[],@mean)];
std_changedpix18 = [ud, accumarray(iy,planform_overlap_ob_18(:,2),[],@std)];
%Phi
mean_Phi18 = [ud, accumarray(iy,planform_overlap_ob_18(:,3),[],@mean)];
std_Phi18 = [ud, accumarray(iy,planform_overlap_ob_18(:,3),[],@std)];
%Corr
mean_corr18 = [ud, accumarray(iy,planform_overlap_ob_18(:,4),[],@mean)];
std_corr18 = [ud, accumarray(iy,planform_overlap_ob_18(:,4),[],@std)];

%Treatment Data
[ud,ix,iy] = unique(planform_overlap_ob_19(:,1));
%changed pix
mean_changedpix19 = [ud, accumarray(iy,planform_overlap_ob_19(:,2),[],@nanmean)];
std_changedpix19 = [ud, accumarray(iy,planform_overlap_ob_19(:,2),[],@nanstd)];
%Phi
mean_Phi19 = [ud, accumarray(iy,planform_overlap_ob_19(:,3),[],@nanmean)];
std_Phi19 = [ud, accumarray(iy,planform_overlap_ob_19(:,3),[],@nanstd)];
%Corr
mean_corr19 = [ud, accumarray(iy,planform_overlap_ob_19(:,4),[],@nanmean)];
std_corr19 = [ud, accumarray(iy,planform_overlap_ob_19(:,4),[],@nanstd)];

%Now lets plot the data

%arrays for plotting
%changed pix
array18_changedpix = [(1:560); mean_changedpix18(:,2)'; std_changedpix18(:,2)'];
cols = any(isnan(array18_changedpix),1);
array18_changedpix(:,cols) = [];

array19_changedpix = [(1:560);mean_changedpix19(:,2)'; std_changedpix19(:,2)'];
cols = any(isnan(array19_changedpix),1);
array19_changedpix(:,cols) = [];

%Phi
array18_phi = [(1:560); mean_Phi18(:,2)'; std_Phi18(:,2)'];
cols = any(isnan(array18_phi),1);
array18_phi(:,cols) = [];

array19_phi = [(1:560);mean_Phi19(:,2)'; std_Phi19(:,2)'];
cols = any(isnan(array19_phi),1);
array19_phi(:,cols) = [];

%Corr
array18_corr = [(1:560); mean_corr18(:,2)'; std_corr18(:,2)'];
cols = any(isnan(array18_corr),1);
array18_corr(:,cols) = [];

array19_corr = [(1:560);mean_corr19(:,2)'; std_corr19(:,2)'];
cols = any(isnan(array19_corr),1);
array19_corr(:,cols) = [];

%fill standard deviation
%changed pix
y18_chpix = array18_changedpix(2,:); % your mean vector;
x18_chpix = array18_changedpix(1,:);
std18_chpix = array18_changedpix(3,:);
curve1_18_chpix = y18_chpix + std18_chpix;
curve2_18_chpix = y18_chpix - std18_chpix;

y19_chpix = array19_changedpix(2,:); % your mean vector;
x19_chpix = array19_changedpix(1,:);
std19_chpix = array19_changedpix(3,:);
curve1_19_chpix = y19_chpix + std19_chpix;
curve2_19_chpix = y19_chpix - std19_chpix;

%Phi
y18_phi = array18_phi(2,:); % your mean vector;
x18_phi = array18_phi(1,:);
std18_phi = array18_phi(3,:);
curve1_18_phi = y18_phi + std18_phi;
curve2_18_phi = y18_phi - std18_phi;

y19_phi = array19_phi(2,:); % your mean vector;
x19_phi = array19_phi(1,:);
std19_phi = array19_phi(3,:);
curve1_19_phi = y19_phi + std19_phi;
curve2_19_phi = y19_phi - std19_phi;

%corr
y18_corr = array18_corr(2,:); % your mean vector;
x18_corr = array18_corr(1,:);
std18_corr = array18_corr(3,:);
curve1_18_corr = y18_corr + std18_corr;
curve2_18_corr = y18_corr - std18_corr;

y19_corr = array19_corr(2,:); % your mean vector;
x19_corr = array19_corr(1,:);
std19_corr = array19_corr(3,:);
curve1_19_corr = y19_corr + std19_corr;
curve2_19_corr = y19_corr - std19_corr;

%%Plot the data
fig = figure()
%plot changed pixels
subplot(3,1,1)
plot(x18_chpix, y18_chpix, 'b', 'LineWidth', 2)
hold on
plot(x19_chpix, y19_chpix, 'g', 'LineWidth', 2)
patch([x18_chpix fliplr(x18_chpix)], [curve1_18_chpix fliplr(curve2_18_chpix)], 'b')
patch([x19_chpix fliplr(x19_chpix)], [curve1_19_chpix fliplr(curve2_19_chpix)], 'g')
alpha(0.15)
plot(x18_chpix, y18_chpix, 'b', 'LineWidth', 2)
plot(x19_chpix, y19_chpix, 'g', 'LineWidth', 2)
%ylim([0 1])
%xlim([0 300])
grid on
grid minor
ylabel('number of changed pixels')
xlabel('measurement window (hr)') 
legend('control mean', 'treatment mean', 'control stdev', 'treatment stdev')
%plot Phi 
subplot(3,1,2)
plot(x18_phi, y18_phi, 'b', 'LineWidth', 2)
hold on
plot(x19_phi, y19_phi, 'g', 'LineWidth', 2)
patch([x18_phi fliplr(x18_phi)], [curve1_18_phi fliplr(curve2_18_phi)], 'b')
patch([x19_phi fliplr(x19_phi)], [curve1_19_phi fliplr(curve2_19_phi)], 'g')
alpha(0.15)
plot(x18_phi, y18_phi, 'b', 'LineWidth', 2)
plot(x19_phi, y19_phi, 'g', 'LineWidth', 2)
%ylim([0 1])
%xlim([0 100])
grid on
grid minor
ylabel('Phi (-)')
xlabel('measurment window (hr)') 
legend('control mean', 'treatment mean', 'control stdev', 'treatment stdev')
%plot correlation
subplot(3,1,3)
plot(x18_corr, y18_corr, 'b', 'LineWidth', 2)
hold on
plot(x19_corr, y19_corr, 'g', 'LineWidth', 2)
patch([x18_corr fliplr(x18_corr)], [curve1_18_corr fliplr(curve2_18_corr)], 'b')
patch([x19_corr fliplr(x19_corr)], [curve1_19_corr fliplr(curve2_19_corr)], 'g')
alpha(0.15)
plot(x18_corr, y18_corr, 'b', 'LineWidth', 2)
plot(x19_corr, y19_corr, 'g', 'LineWidth', 2)
%ylim([0 1])
%xlim([0 100])
grid on
grid minor
ylabel('correlation (-)')
xlabel('measurment window (hr)') 
legend('control mean', 'treatment mean', 'control stdev', 'treatment stdev')
%patch([min(xlim) max(xlim) max(xlim) min(xlim)], [ybars(1) ybars(1), ybars(2) ybars(2)], [0.8 0.8 0.8])
set(gcf, 'PaperUnits', 'inches');
y_width=7.25;x_width=9.125
set(gcf, 'PaperPosition', [0 0 x_width y_width]);
saveas(fig, 'allflow_planform_change.pdf')
