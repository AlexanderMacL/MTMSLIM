% Example using SLIMwrap function to analyse a single image and show results

clear all
close all

% [imagenames{1},imagepath] = uigetfile({'*.bmp','Bitmap images (*.bmp)'},'SLIM image');
% [clnm,clpth] = uigetfile({'*.txt','SLIM calibration files (*.txt)'},'SLIM calibration file');
% calpath = [clpth, clnm];

% More detailed explanations of these parameters are given in SLIMwrap.m
imagepath = '';
imagenames{1} = 'ExampleImage.bmp';
calpath = 'ExampleCalibration.txt';
maxdistance = 0; % values of F whose corresponding D is greater than maxdistance are set to zero, can be used to suppress noise
minfilm = 0; % if nonzero, truncates the cal curve at this lower bound before matching
maxfilm = 250; % if nonzero, truncates the cal curve at this upper bound before matching
masksize = 60; % second-nearest colour is determined by masking this number of nanometres of the calibration curve either side of the nearest colour, and so on.
calcurverefit = true; % refit the calibration curve based on pixel values in the image
calcurveoffsetprescale = 120; % fudge factor for the shape of the cal curve, between 0 and 100 is usually OK, cannot be less than -1* the minimum value of the cal curve on any channel, larger value implies greater skew
numorders = 2; % number of 'nearest colours' e.g. 2 considers both nearest and second-nearest colours when minimising 2nd derivative
numaltmaps = 2; % number of alternative maps based on the assumption that the first 5 pixels correspond to the nearest, 2nd-nearest-colour film thickness etc. Cannot be greater than numorders
plotfigures = true; % set this to false to stop seeing the cal curve refit plots
tangentcorrection = false; % experimental feature, usually causes bad things to happen, leave false
circlefitting = true; % if true, finds the contact circle, limits the 
    % field of view to a rectangle containing it, and sets values of F
    % outside the circle to zero. This may increase analysis time slightly
avgblocksize = 3;
sampcircrad = 0;
sampcirccent = [768,1024]./2;

[A,AA,F,D,R,Favg] = SLIMwrap(imagepath,imagenames,calpath,maxdistance,minfilm,maxfilm,circlefitting,avgblocksize,sampcircrad,sampcirccent,calcurverefit,calcurveoffsetprescale,numorders,numaltmaps,masksize,tangentcorrection,plotfigures);

mxx = max(cellfun(@(x) max(x,[],'all'),F));
mnn = min(cellfun(@(x) min(x,[],'all'),F));
mxd = max(cellfun(@(x) max(x,[],'all'),D));
mnd = min(cellfun(@(x) min(x,[],'all'),D));
% for i=1:length(F)
%     writematrix(F{i},[fnames{i}(1:end-4),'_FILMTHICKNESS.csv']);
%     writematrix(D{i},[fnames{i}(1:end-4),'_DISTANCE.csv']);
% end

figure();
montage(D);
colorbar;
caxis([mnd, mxd]);
title("Distance")
saveas(gcf,[imagepath 'Distance_wrap.png']);

figure();
montage(A);
title("Original Image")
saveas(gcf,[imagepath 'Mapper_wrap.png']);

figure();
montage(AA);
title("Cropped and Averaged Image")
saveas(gcf,[imagepath 'MapperCropped_wrap.png']);

figure();
montage(F);
set(gca,'FontSize',28)
colormap(gca,turbo());
c = colorbar;
c.Label.String = 'Film Thickness [nm]';
caxis([mnn, mxx]);
title("Film Thickness")
saveas(gcf,[imagepath 'Film_wrap.png']);

figure();
montage(R);
title("Reconstructed Colour Image")
saveas(gcf,[imagepath 'Reconstructed_wrap.png']);


