% Example using MTMSLIM function

clear all
close all

% Use the UI file picker to find your mtmd file, which allows the SLIM
% images to be located assuming they are in a folder in the same
% directory with the same name appended by 'slim images'
[testname,testpath] = uigetfile({'*.mtmd;*.etmd','MTM/ETM data files (*.mtmd;*.etmd)'});

% Alternatively, specify the path and filename yourself
% testpath = 'Tests\Example\'; % directory containing mtmd file
% testname = 'Example.mtmd'; % name of mtmd/etmd file including extension

% More detailed explanations of these parameters are given in SLIMwrap.m
maxdistance = 0; % values of F whose corresponding D is greater than maxdistance are set to zero, can be used to suppress noise
minfilm = -20; % if nonzero, truncates the cal curve at this lower bound before matching
maxfilm = 200; % if nonzero, truncates the cal curve at this upper bound before matching
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
avgblocksize = 5;
sampcircrad = 0;
sampcirccent = [768,1024]./2;

[A,AA,F,D,R,fnames,Favg,mapi,J] = MTMSLIM(testpath,testname,maxdistance,minfilm,maxfilm,circlefitting,avgblocksize,sampcircrad,sampcirccent,calcurverefit,calcurveoffsetprescale,numorders,numaltmaps,masksize,tangentcorrection,plotfigures);

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
set(gca,'FontSize',48);
c = colorbar;
caxis([mnd, mxd]);
c.Label.String = 'Distance (L2) [-]';
% title("Distance")
set(gcf,'Position',[540,200,1200,760]);
saveas(gcf,[testpath testname(1:end-5) '_Distance.png']);

figure();
montage(A);
title("Original Image")
saveas(gcf,[testpath testname(1:end-5) '_Mapper.png']);

figure();
montage(AA);
title("Cropped and Averaged Image")
saveas(gcf,[testpath testname(1:end-5) '_MapperCropped.png']);

figure();
montage(F);
set(gca,'FontSize',48);
colormap(gca,turbo());
c = colorbar;
caxis([mnn, mxx]);
caxis([-40,240]);
c.Label.String = 'Film Thickness [nm]';
% title("Film Thickness")
set(gcf,'Position',[540,200,1200,760]);
saveas(gcf,[testpath testname(1:end-5) '_Film.png']);

figure();
montage(R);
title("Reconstructed Colour Image")
saveas(gcf,[testpath testname(1:end-5) '_Reconstructed.png']);

figure();
plot(Favg,'-o');
xlabel('Mapper image # [-]')
ylabel('Mean Film Thickness (nm)')
saveas(gcf,[testpath testname(1:end-5) '_FilmAvg.png']);