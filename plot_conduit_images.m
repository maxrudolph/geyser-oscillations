clear;
close all;

path = '/Volumes/GeyserData/NSFGeyserProject/HighSpeed/11-22-2024/F7/';
framerate = 500;

max_frame = 1000;
% compute mean over first X frames
for frame=1:max_frame
    filename = sprintf([path 'img_%06d.png'],frame);
    im = imread(filename);
    if frame == 1
        immean = double(im)/max_frame;
    else
        immean = immean + double(im)/max_frame;
    end
end
%%
immean_u8 = uint8(immean);
figure, imshow(immean_u8);

%% read images for plotting

i=1;
im1 = {};
output_times = [];
for frame=4000:1000:12000
    filename = sprintf([path 'img_%06d.png'],frame);
    im1{i} = imread(filename);
    %% add black bars
    im1{i}(1:5,:) = 0;
    im1{i}(end-4:end,:) = 0;

    if i == 1        
        imcat = im1{1};
    else        
        imcat = cat(1,imcat,im1{i});
    end
    output_times(i) = frame/framerate;
    i=i+1;
    
end
output_times = output_times - output_times(1);
ntot = i-1;

immean_bw = rgb2gray(immean);
% combine images
% test = cat(1,im)
imshow(imcat);
imbw = rgb2gray(imcat);
enhanced = adapthisteq(imbw);
figure,
imshow(flipud(enhanced'));

for i=1:ntot
    text(1/ntot*(i-1)+.03,0.05,sprintf('%.0f s',output_times(i)),'Fontsize',24,'Units','normalized','Color','w')
end

exportgraphics(gcf,'bottom_conduit_series.pdf');