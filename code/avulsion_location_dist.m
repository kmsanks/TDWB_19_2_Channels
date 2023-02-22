clear all; close all

cd '../data'
control = csvread('control_avulsions.csv');
treatment = csvread('treatmeant_avulsions.csv');

%apex 
xentrance_18 = 109;%x grid node location of the entrance channel
yentrance_18 = 271;%y grid node location of the entrance channel

xentrance_19 = 214;%x grid node location of the entrance channel (x is down dip)
yentrance_19 = 397;%y grid node location of the entrance channel (y is strike)

La18 = sqrt(((control(:,1)-xentrance_18)*5).^2 + ((control(:,2)-yentrance_18)*5).^2);
La19 = sqrt(((treatment(:,1)-xentrance_19)*5).^2 + ((treatment(:,2)-yentrance_19)*5).^2);

figure()
histogram(La18, 'FaceColor','blue','BinWidth', 100)
hold on
histogram(La19, 'FaceColor','green', 'FaceAlpha',0.5, 'BinWidth',100)

mean(La18)
mean(La19)
