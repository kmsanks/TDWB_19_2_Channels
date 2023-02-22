%%lets calculate slope above and below backwater point
%we want to remove first 100 mm or so from the entrance channel because the
%elevations are often lower there
slope_dat18 = backwater_mat_18(10:end, :, 1:2:560);
slope_dat19 = backwater_mat_19(10:end, :, :);

%slope
slope_array18 = [];
for i = 1:size(slope_dat18, 2)
    dat = slope_dat18(:,:,i);
    if isnan(backwater_dist18(i))
        slope_array18(:,i) = [NaN;NaN];
    else
        idx1 = 1;
        idx2 = find(dat(:,1) == backwater_dist18(i));
        idx3 = idx2 + 1;
        idx4 = find(dat(:,1) > backwater_dist18(i));
        idx4 = idx4(end);
        fx = fitlm(dat(idx1:idx2, 1), dat(idx1:idx2,2));  %fit linear model for beginning of data 
        slope_begin = table2array(fx.Coefficients(2,1));
        fx = fitlm(dat(idx3:idx4, 1), dat(idx3:idx4,2));   %fit linear model for end of data
        slope_end = table2array(fx.Coefficients(2,1));
        slope_array18(:,i) = [slope_begin,slope_end]';    %save beginning slope and end slope
    end
end

slope_array19 = [];
for i = 1:size(slope_dat19, 2)
    dat = slope_dat19(:,:,i);
    if(isnan(dat(:,2)))
         slope_array19(:,i) = [NaN,NaN]';
    else
        idx1 = 1;
        idx2 = find(dat(:,1) == backwater_dist19(i));
        idx3 = idx2 + 1;
        idx4 = find(dat(:,1) > backwater_dist19(i));
        idx4 = idx4(end);
        fx = fitlm(dat(idx1:idx2, 1), dat(idx1:idx2,2));  %fit linear model for beginning of data 
        slope_begin = table2array(fx.Coefficients(2,1));
        fx = fitlm(dat(idx3:idx4, 1), dat(idx3:idx4,2));   %fit linear model for end of data
        slope_end = table2array(fx.Coefficients(2,1));
        slope_array19(:,i) = [slope_begin,slope_end]';    %save beginning slope and end slope
    end
end

%%Here is a method to find slope breaks in the channels
% %we want to remove first 100 mm or so from the entrance channel because the
% %elevations are often lower there
% slope_dat18 = backwater_mat_18(10:end, :, 1:2:560);
% slope_dat19 = backwater_mat_19(10:end, :, :);
% %find the slope break point in the channels
% breakpts_18 = [];
% for i = 1:size(slope_dat18,3)
%     sbreak = slmengine(slope_dat18(:,1,i),slope_dat18(:,2,i),'degree',1,'plot','on','knots',3,'interiorknots','free'); %calculates slope break
%     breakpts_18(:,i) = sbreak.knots;
% end
% 
% 
% breakpts_19 = [];
% for i = 1:size(slope_dat19,3)
%     if isnan(slope_dat19(:,2,i))
%         breakpts_19(:,i) = [NaN,NaN,NaN]';
%     else
%         sbreak = slmengine(slope_dat19(:,1,i),slope_dat19(:,2,i),'degree',1,'plot','on','knots',3,'interiorknots','free'); %calculates slope break
%         breakpts_19(:,i) = sbreak.knots;
%     end
% end


% %lets calculate slope on the break points
% slope_array18 = [];
% for i = 1:size(breakpts_18, 2)
%     dat = slope_dat18(:,:,i);
%     idx1 = find(dat(:,1) == breakpts_18(1,i));
%     idx2 = find(dat(:,1) >= breakpts_18(2,i)); %the slope break is not necessarily at a data point, 
%     idx2 = idx2(1);                            %so we need the first point greater than this location
%     idx3 = find(dat(:,1) == breakpts_18(3,i));
%     fx = fitlm(dat(idx1:idx2, 1), dat(idx1:idx2,2));  %fit linear model for beginning of data 
%     slope_begin = table2array(fx.Coefficients(2,1));
%     fx = fitlm(dat(idx2:idx3, 1), dat(idx2:idx3,2));   %fit linear model for end of data
%     slope_end = table2array(fx.Coefficients(2,1));
%     slope_array18(:,i) = [slope_begin,slope_end]';    %save beginning slope and end slope
% end
% 
% 
% %lets calculate slope on the break points
% slope_array19 = [];
% for i = 1:size(breakpts_19, 2)
%     dat = slope_dat19(:,:,i);
%     if(isnan(dat(:,2)))
%          slope_array18(:,i) = [NaN,NaN]';
%     else
%         idx1 = find(dat(:,1) == breakpts_19(1,i));
%         idx2 = find(dat(:,1) >= breakpts_19(2,i)); %the slope break is not necessarily at a data point, 
%         idx2 = idx2(1);                            %so we need the first point greater than this location
%         idx3 = find(dat(:,1) == breakpts_19(3,i));
%         fx = fitlm(dat(idx1:idx2, 1), dat(idx1:idx2,2));  %fit linear model for beginning of data 
%         slope_begin = table2array(fx.Coefficients(2,1));
%         fx = fitlm(dat(idx2:idx3, 1), dat(idx2:idx3,2));   %fit linear model for end of data
%         slope_end = table2array(fx.Coefficients(2,1));
%         slope_array19(:,i) = [slope_begin,slope_end]';    %save beginning slope and end slope
%     end
% end
% 
% figure
% plot(1:2:560,slope_array18(1,:), 'bo')
% hold on
% plot(1:2:560,slope_array19(1,:), 'go')
% plot(1:2:560,slope_array18(2,:), 'bx')
% plot(1:2:560,slope_array19(2,:), 'gx')
% 
% nanmedian(slope_array18(1,:))
% nanmedian(slope_array19(1,:))
% nanmedian(slope_array18(2,:))
% nanmedian(slope_array19(2,:))
% 
% figure()
% x = [ones(size(slope_array18(1,:))); 2*ones(size(slope_array18(2,:))); 3*ones(size(slope_array19(1,:))); 4*ones(size(slope_array19(2,:)))];
% g = [slope_array18(1,:); slope_array18(2,:); slope_array19(1,:); slope_array19(2,:)];
% boxplot(x,g) %, 'labels', {'control start', 'control_end', 'treatment start', 'treatment end'});