%% 
% Code information
% Code to Calculate sigma ss for range of measurment windows provided by
% data set(ts) with delta time between measurements(dt) and an accomodation production rate (drift).
% ts matrix should be aligned as Space,Time
function [t,sss,TC,kappa_st,kappa_lt] = sigma_ss(ts,dt,drift);
%% Definitions
nc=size(ts,1);%Number of spatial locations
nr=size(ts,2);%Number of time steps
%drift=0.25;%Mean drift per time step (mm/hr)
sss=[];%sigma ss values
sss_5 = [];
sss_95 = [];
%dt = 1;%delta t of time series (hr)
t = [dt:dt:dt*(nr-1)];
%% Main loop for sigma ss calculations
for i=1:1:(nr-1)%loop with changing window size
    Utspan=[];%ratios of sedimentation/subsidence
    for j=1:nr
        line1=ts(:,j);%topography at time t=j
        if (i+j)>nr
            break
        end
        line2=ts(:,(i+j));%topography at time t=j+i
        l1m2=line2-line1;%difference in line 2 and line 1 (i.e. sedimentation)
        bslevelrise=i*drift;%mean drift for this time period
        U=l1m2./bslevelrise;%ratio of sedimentation/subsidence
        Utspan=[Utspan;std(U, 'omitnan')];
    end
    sttspan = mean(Utspan, 'omitnan');
    sttspan5 = prctile(Utspan,5);
    sttspan95 = prctile(Utspan,95);
    sss=[sss;sttspan];
    sss_5 = [sss_5;sttspan5];
    sss_95 = [sss_95;sttspan95];
end
%% crop data
t = t(1:floor(max(size(t))/2));
sss = sss(1:floor(max(size(sss))/2));
sss_5 = sss_5(1:floor(max(size(sss_5))/2));
sss_95 = sss_95(1:floor(max(size(sss_95))/2));
sss_95eb = sss_95 - sss;
sss_5eb = sss - sss_5;
%% compensation index, kappa_st, kappa_lt calcs
lnt = log(t);
lnsss = log(sss);
dp = 0.001;
lnt_int = min(lnt):dp:max(lnt);
lnsss_int = interp1(lnt,lnsss,lnt_int);
FC = findchangepts(lnsss_int,'Statistic','linear','MaxNumChanges',1);
C1 = [exp(lnt_int(FC(1)));exp(lnt_int(FC(1)))];
YC = [min(sss_5) max(sss_95)];
TC = C1(1)
[c index] = min(abs(t-TC(1)));
%Fit = ezfit(t(1:(index-1)),sss(1:(index-1)),'power;log');
Fit = ezfit(t(1:round(index/2)),sss(1:round(index/2)),'power;log');
co = Fit.m;
kappa_st = -1*co(2)%short-term compensation index, should be between 0 and 1
Fit2 = ezfit(t((index+1):end),sss((index+1):end),'power;log');
co = Fit2.m;
kappa_lt = -1*co(2)%long-term compensation index, should be ~1
%% Figures
figure(1)
errorbar(t,sss,sss_5eb,sss_95eb,'.k')
hold on
loglog(t,sss,'md')
loglog(C1,YC,'r')
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
ylabel('sigma ss')
xlabel('time window (hr)')