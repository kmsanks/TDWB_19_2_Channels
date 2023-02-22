clear all

cd '..\data'
load('CM_18.mat')
load('CM_19.mat')
cd '.\images18'
files = dir();
%set(gca, 'YDir', 'reverse')
figure()
for i = 3:length(files)
    img = imread(files(i).name);
    subplot(1,2,1)
    imagesc(img)
    title(num2str(i-2))
    img2 = imread(files(i+3).name);   
    subplot(1,2,2)
    imagesc(img2)
    title(num2str(i+3-2))
    pause;
end

cd '..\images19'
files = dir();
%set(gca, 'YDir', 'reverse')
figure()
for i = 3:length(files)
    img = imread(files(i).name);
    imagesc(img)
    pause;
end
