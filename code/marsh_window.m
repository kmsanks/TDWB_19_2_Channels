%% Calculate distance to the marsh window
% determine elevations outside channel
FFM_18 = abs(CM_18-1);
FFM_19 = abs(CM_19-1);

%if no channel map in treatment
for i = 1:size(FFM_19,3)
    if sum(sum(FFM_19(:,:,i))) == 560250;
        FFM_19(:,:,i) = NaN;
    end
end

% Remove area inside channels 
% control
zff18 = [];
for i =1:size(FFM_18,3)
  zff=(FFM_18(:,:,i)).*(z18(:,:,i));
  zff(zff == 0) = NaN; %0 is inside channel
  zff18(:,:,i) = zff; % matrix of elevations outside channel
end
% treatment
zff19 = NaN(size(FFM_19));
for i = 1:size(FFM_19,3)
  % if no channel map, matrix should be NaN
    tmp = FFM_19(:,:,i);
    if  isnan(sum(tmp(:)))
        zff19(:,:,i) = tmp;
    else % make matrix of elevations outside channel 
        zff=(FFM_19(:,:,i)).*(z19(:,:,i));
        zff(zff == 0) = NaN; %0 is inside channel
        zff19(:,:,i) = zff;
    end
end

% loop to calculate distance to marsh window
marsh_elev18 = []; %initialize matrix
marsh_dist18 = []; %initialize matrix
dist = 0:10:3100; %radial loop distances
marsh_mat_18 = NaN(length(dist)-1, 2, size(FFM_18,3));
for i = 1:size(FFM_18,3)
    marsh_elev = []; %initialize matrix
    marsh_dist = []; %initialize matrix
    for k = 1:(length(dist)-1) 
        idx = dd18(dd18 >= dist(k) & dd18 < dist(k+1));
        radial_dd = dd18 >= dist(k) & dd18 < dist(k+1);
        is_shot = FFM_18(:,:,i).*radial_dd;
        is_shot(is_shot == 0.) = NaN;
        z = zff18(:,:,i).*is_shot; %elevation relative to sea level
        %now we want to determine if at least 50% of the elevations are in
        %the marsh window
        el = prctile(z(:),50); %median elevation of transect
        if el <= 5 %if it is in the marsh window, save elevation and distance
           marsh_elev = [marsh_elev; el];
           marsh_dist = [marsh_dist; dist(k)];
        end
        marsh_mat_18(k,1,i) = dist(k); %save all data to matrix
        marsh_mat_18(k,2,i) = el; % save all data to matrix
    end
    if isempty(marsh_elev)
       marsh_elev18(i) = NaN; %never makes it to marsh window
       marsh_dist18(i) = NaN; %never makes it to marsh window
    else
       marsh_elev18(i) = marsh_elev(1); %first transect in marsh window
       marsh_dist18(i) = marsh_dist(1); %first transect in marsh window
    end 
end

% treatment 
marsh_elev19 = [];
marsh_dist19 = [];
marsh_mat_19 = NaN(length(dist)-1, 2, size(FFM_19,3));
for i = 1:size(FFM_19,3)
    marsh_elev = [];
    marsh_dist = [];
    for k = 1:(length(dist)-1) 
        idx = dd19(dd19 >= dist(k) & dd19 < dist(k+1));
        radial_dd = dd19 >= dist(k) & dd19 < dist(k+1);
        is_shot = FFM_19(:,:,i).*radial_dd;
        [row col] = find(~isnan(is_shot));
        if isempty(row)
            marsh_elev19(i) = NaN;
            marsh_dist19(i) = NaN;
            marsh_mat_19(k,1,i) = NaN;
            marsh_mat_19(k,2,i) = NaN;
        else
            is_shot(is_shot == 0.) = NaN;
            z = zff19(:,:,i).*is_shot; %elevation relative to sea level
            %now we want to determine if at least 50% of the elevations are in
            el = prctile(z(:),50); %median elevation of transect
            if el <= 5 %if it is in the marsh window, save elevation and distance
               marsh_elev = [marsh_elev; el];
               marsh_dist = [marsh_dist; dist(k)];
            end
            marsh_mat_19(k,1,i) = dist(k); %save to matrix
            marsh_mat_19(k,2,i) = el; %save to matrix
        end
        if isempty(marsh_elev)
           marsh_elev19(i) = NaN; %never makes it to marsh window
           marsh_dist19(i) = NaN; %never makes it to marsh window
        else
           marsh_elev19(i) = marsh_elev(1); %first transect in marsh window
           marsh_dist19(i) = marsh_dist(1); %first transect in marsh window
        end 
    end
end


figure
plot(backwater_point18/1000, marsh_dist18(1:559)/1000., 'bo')
hold on
%plot(backwater_dist19/1000, marsh_ent19_50/1000., 'go')
plot(mdl_treatmarsh, 'color', 'green', 'marker', 'o')
plot(0:1:3,0:1:3, 'k-')
xlabel('distance from apex to backwater point (m)')
ylabel('distance from apex to marsh entrance (m)')
xlim([0 2])
ylim([0 2])
legend('control', 'treatment', 'data', 'linear fit', 'CI', '1:1 line')

fig2 = figure();
subplot(3,1,1)
plot(1:2:560, backwater_dist18(1:2:560)/1000, 'bo')
hold on
plot(hour18(1:2:length(hour18)), nan_bw18(1:2:length(hour18)), 'bx')
plot(1:2:560, backwater_dist19/1000, 'go')
plot(hour19', nan_bw19, 'gx')
%plot(hour19*2, nan_bw19, 'gx') %these are the same time steps as the ones
%with no channel maps
plot(hour_nomaps19*2, bw_nomaps19, 'kx')
xlabel('time (hours)')
ylabel('distance from entrance to backwater location (m)')
xlim([0 560])
legend('control', 'no bw control', 'treatment', 'no bw treatment', 'no chan map treatment')
subplot(3,1,2)
plot(1:2:560, backwater_length18(1:2:560)/1000, 'bo')
hold on
plot(1:2:560, backwater_length19/1000, 'go')
plot(hour_nomaps19*2, bw_nomaps19, 'kx')
plot(hour_nomaps19*2, bw_nomaps19, 'kx')
xlabel('time (hours)')
ylabel('backwater length (m)')
xlim([0 560])
legend('control', 'treatment', 'no chan map treatment')
subplot(3,1,3)
plot(backwater_dist18(1:2:560)/1000, marsh_ent18_50/1000., 'bo')
hold on
%plot(backwater_dist19/1000, marsh_ent19_50/1000., 'go')
plot(mdl_treatmarsh, 'color', 'green', 'marker', 'o')
plot(0:1:3,0:1:3, 'k-')
xlabel('distance from entrance to backwater location (m)')
ylabel('distance from entrance to radial marsh transect (m)')
xlim([0 2])
ylim([0 2])
legend('control', 'treatment', 'data', 'linear fit', 'CI', '1:1 line')
% subplot(3,1,3)
% plot(backwater_dist18(1:2:560)/1000, marsh_exit18/1000., 'bo')
% hold on
% plot(backwater_dist19/1000, marsh_exit19/1000., 'go')
% plot(0:1:3,0:1:3, 'k-')
% xlabel('distance from entrance to backwater location (m)')
% ylabel('distance from entrance to furthest aerial pixel (m)')
% xlim([0 2])
% %ylim([0 1])
% legend('control', 'treatment', '1:1 line')
y_width=7.25;x_width=9.125
set(gcf, 'PaperPosition', [0 0 x_width y_width]);
saveas(fig2, 'backwater_stats_v2.pdf')

