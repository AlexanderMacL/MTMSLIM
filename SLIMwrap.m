% MTM SLIM Calculation abstracted to not require a profile
% (rig-independent)

% Author: Alexander MacLaren
% Last updated: 03/04/2024

% Premise:
% Read in images specified
% Find circles, crop, average
% Calcurverefit
% Apply expandingbox and get pixel sequence
% For first pixel in sequence, find nearest film, and next-nearest etc up to
% numorders
% Colour match for each pixel in sequence, for each branch in numorders,
% each time updating the cal curve margins in each branch
% Store all the things for every branch, and at the end choose the
% minimum-distance branch to go with


% Arguments:
%   imagepath is a string or character vector containing the directory in
%       which the images are located, including trailing /
%   imagenames is a cell array of strings or character vectors containing
%       the names of the images in imagepath
%   calpath is the full path of the calibration file
%   maxdistance is a number, values of F whose corresponding D is greater than
%       maxdistance are set to zero, can be used to suppress noise
%   minfilm (optional) is a number, values of F less than minfilm are
%       disallowed and the next-nearest film thickness is found, can be
%       used to suppress noise - if minfilm is missing or zero, no minimum
%       is imposed
%   maxfilm (optional) is a number, values of F greater than maxfilm are
%       disallowed and the next-nearest film thickness is found, can be
%       used to suppress noise - if maxfilm is missing or zero, no maximum
%       is imposed
%   distancemetric (optional) is a string or character vector specifying the distance
%       metric, supported values are 'L1' (linear) and 'L2' (least squares)
%       if curverefitting is enabled, distancemetric is ignored and L2 is
%       used
%   circlefitting (optional) is a boolean. If true, finds the contact
%       circle, limits the field of view to a rectangle containing it, and
%       sets values of F outside the circle to zero. This may increase
%       analysis time slightly (defaults to false). Averaging is only
%       conducted when circlefitting is false.
%   avgblocksize (optional) is an odd positive integer specifying the side
%       length in pixels of the averaging grid to be tiled. A value of 1
%       results in no averaging. (defaults to 1)
%   sampcircrad (optional) is an even positive integer specifing the radius
%       of the sampling circle, and must be greater than avgblocksize.
%       If zero, circle cropping is not performed. (defaults to 384) 
%   sampcirccent (optional) is a 1x2 row vector of co-ordinates of the
%       centre of the sampling circle (defaults to the centre of the image)
%   calcurverefit (optional) is a boolean, which if true prompts a refit
%       of the calibration curve to the image, which is useful if the
%       saturation of the image doesn't match the calibration, although
%       this significantly increases analysis time
%   calcurveoffsetprescale (optional) is a positive real number which
%       parameterises the skew in the calibration curve and may be tweaked
%       to improve the quality of the calibration curve fit.
%       Defaults to 120
%   numorders (optional) is an integer specifying the number of successive
%       next-nearest branches of the calibration curve to consider - a
%       value of 1 implies strict nearest-colour, defaults to 4
%   numaltmaps (optional) is an integer specifying the number of
%       simultaneous film thickness maps constructed based on the
%       next-nearest fit to the first 5 pixels in the expandingbox
%       sequence. This enables the fit to consider alternatives if the
%       nearest colour film thickness to the first 5 pixels is incorrect.
%       The minimum-total-distance map is the one finally selected.
%       Defaults to 1, cannot be greater than numorders
%   masksize (optional) is an integer specifying the half-width in
%       nanometres of the mask which is applied to the calibration curve to
%       force the next-nearest colour fit - defaults to 40, approx. a
%       quarter blue-intensity-curve-wavelength
%   tangentcorrection (optional) is a boolean which if true switches from
%       traverse-perpendicular second-derivative minimisation to
%       traverse-tangential, in cases where the traverse is not interrupted
%       by an edge (defaults to false)
%   plotfigures (optional) is a boolean, and enables the plotting of
%      various 3D plots to show the relationship between the calibration
%      curve, the fit, and the data, at different points in the analysis

% Return values:
%   A is a cell array containing RGB images showing contact
%   AA is a cell array containing cropped RGB images
%   F is a cell array containing contact film thickness maps
%   D is a cell array containing distance maps
%   R is a cell array containing reconstructed images
%   fnames is a cell array of strings containing the image filenames
%   Favg is a vector containing the average film thickness


function [A,AA,F,D,R,Favg] = SLIMwrap(imagepath, imagenames, calpath, ...
    maxdistance, minfilm, maxfilm, circlefitting, avgblocksize, ...
    sampcircrad, sampcirccent, calcurverefit, calcurveoffsetprescale, ...
    numorders, numaltmaps, masksize, tangentcorrection, plotfigures)

CIRCLEFIT_METHOD = 3;
% 1 : Hough Transform method using imfindcircles
% 2 : Convolution circle fitting using SLIMConvFindCircle
% 3 : Separated-axis circle fitting using FastFeatureFind

if nargin<17
    plotfigures = false;
    if nargin<16
        tangentcorrection = false;
        if nargin<15
            masksize = 40; % 1/4 of a blue wavelength
            if nargin<14
                numaltmaps = 1;
                if nargin<13
                    numorders = 4;
                    if nargin<12
                        calcurveoffsetprescale = 120;
                        if nargin<11
                            calcurverefit = true;
                            if nargin<10
                                sampcirccent = [768,1024]./2;
                                if nargin<9
                                    sampcircrad = 768/2;
                                    if nargin<8
                                        avgblocksize = 1;
                                        if nargin<7
                                            circlefitting = false;
                                            if nargin<6
                                                maxfilm = 0;
                                                if nargin<5
                                                    maxfilm = 0;
                                                end
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end

if (sampcircrad<=avgblocksize && sampcircrad~=0) || mod(sampcircrad,2)~=0
    error(['sampcircrad (' num2str(sampcircrad) ') must be an even number greater than avgblocksize (' num2str(avgblocksize) ')'])
end

if (numorders<numaltmaps)
    error(['numaltmaps (' num2str(numaltmaps) ') cannot be greater than numorders (' num2str(numorders) ')']);
end

B = [];

% extract calibration
C = readmatrix(calpath,'NumHeaderLines',1,'Delimiter','\t');

if calcurverefit
    C(:,2:4) = C(:,2:4) + calcurveoffsetprescale;
    % regress for cal curve
    betaR = regressexpdecay(C(:,1)',C(:,2)',false,true);
    betaG = regressexpdecay(C(:,1)',C(:,3)',false,true);
    betaB = regressexpdecay(C(:,1)',C(:,4)',false,true);
    if (minfilm ~= 0); minf = minfilm; else minf = min(C(:,1)); end
    if (maxfilm ~= 0); maxf = maxfilm; else maxf = max(C(:,1)); end
    CC(:,1) = minf:1:maxf;
    CC(:,2) = expexpdecsin(betaR,CC(:,1));
    CC(:,3) = expexpdecsin(betaG,CC(:,1));
    CC(:,4) = expexpdecsin(betaB,CC(:,1));
    C(:,2:4) = C(:,2:4) - calcurveoffsetprescale;
    CC(:,2:4) = CC(:,2:4) - calcurveoffsetprescale;
else
    CC = C;
end

fnames = cellfun(@(x) [imagepath x],imagenames,'UniformOutput',false);

imagerange = 1:length(fnames);

A = cellfun(@imread, fnames, 'UniformOutput',false);

disp([num2str(length(fnames)) ' mapper images found'])

[AA,AAA,F,D,R,G,X,Y,XX,f,d,r] = deal(cell(1,length(imagerange)));
Favg = zeros(1,length(imagerange));
for i = 1:length(imagerange)
    disp(['Averaging image ' num2str(i) ' of ' num2str(length(imagerange))]);
    s = size(A{i});
    [X{i},Y{i}] = meshgrid(1:s(2),1:s(1));
    XX{i} = cat(3,X{i},Y{i});
    ringthickness = 5; % pixels
    circlefittingsuccess = circlefitting;
    if circlefitting
        if CIRCLEFIT_METHOD == 1 % HOUGH TRANSFORM CIRCLE FITTING
            AAA{i} = imgaussfilt(sum(double(A{i}),3)./(3*255),4); % imgaussfilt blurs fringes
            [circcent,circrad] = imfindcircles(AAA{i},floor([0.1,0.5].*s(2)),'ObjectPolarity','dark','EdgeThreshold',0.1,'Sensitivity',0.95);
            if length(circrad)<1
                warning(['circlefitting for image ' num2str(imagerange(i)) ' unsuccessful'])
                circlefittingsuccess = false;
                AA{i} = A{i};
            elseif length(circrad)>1
        %         warning(['multiple circles detected - using largest'])
                [circrad(1),ci] = max(circrad);
                circcent(1,:) = circcent(ci,:);
            end
        elseif CIRCLEFIT_METHOD == 2 % CONVOLUTION CIRCLE FITTING
            [circrad,circcent(1,:)] = SLIMConvFindCircle(A{i});
            if circrad==0 && all(circcent(1,:)==0)
                warning(['circlefitting for image ' num2str(imagerange(i)) ' unsuccessful'])
                circlefittingsuccess = false;
                AA{i} = A{i};
            end
        else % FASTFEATUREFIND
            [circcent(1,:),circrad] = FastFeatureFind(A{i});
            if circrad==0 && all(circcent(1,:)==0)
                warning(['circlefitting for image ' num2str(imagerange(i)) ' unsuccessful'])
                circlefittingsuccess = false;
                AA{i} = A{i};
            end
            circcent = flip(circcent,2);
        end
    end
    if circlefittingsuccess
        if sampcircrad==0; sampcircrad = circrad; end
        sampcirccent = circcent(1,:);
    end
    if sampcircrad>0
        avgradextent = floor(sampcircrad/avgblocksize);
    else
        avgradextent = floor(s(2)/2/avgblocksize);
    end
    % THIS ONE HAS BEEN FIXED TO BE SYMMETRICAL
    AVG = zeros(2*avgradextent,2*avgradextent,3);
    for p = -avgradextent:avgradextent-1
        for q = -avgradextent:avgradextent-1
            if sampcircrad>0 && ((((p+0.5)*avgblocksize)^2+((q+0.5)*avgblocksize)^2) > sampcircrad^2)
                AVG(p+avgradextent+1,q+avgradextent+1,:) = 0;
            else
                AVG(p+avgradextent+1,q+avgradextent+1,:) = mean(A{i}(sampcirccent(2)+p*avgblocksize+1:sampcirccent(2)+(p+1)*avgblocksize, sampcirccent(1)+q*avgblocksize+1:sampcirccent(1)+(q+1)*avgblocksize, :),[1,2]);
            end
        end
    end
    if sampcircrad>0
        outofcircindex{i} = sum(AVG,3) == 0;
        ringindex{i} = ((XX{i}(:,:,1)-sampcirccent(1)).^2 + (XX{i}(:,:,2)-sampcirccent(2)).^2) < (sampcircrad + ringthickness)^2;
        ringindex{i} = ringindex{i} & ((XX{i}(:,:,1)-sampcirccent(1)).^2 + (XX{i}(:,:,2)-sampcirccent(2)).^2) > (sampcircrad^2);
    end
    AA{i} = uint8(AVG);
    s = size(AA{i});
    b{i} = double(AA{i});
    if sampcircrad>0; b{i}(repmat(outofcircindex{i},1,1,3)) = NaN; end
    b{i} = reshape(b{i},[s(1)*s(2),3]);
    B = [B;b{i}];
end
if calcurverefit
    [CC]  = CalCurveRefit3MTM(C,CC,betaR,betaG,betaB,B,minfilm,maxfilm,calcurveoffsetprescale,plotfigures,20000);
end
for i = 1:length(imagerange)
    disp(['Analysing image ' num2str(i) ' of ' num2str(length(imagerange))]);
    outofcircindex{i} = sum(AA{i},3) == 0;
    s = size(AA{i},[1,2]);
    [linx,linperp,linintrup] = expandingBox(s,floor(s./2),outofcircindex{i});
    
    if (~tangentcorrection)
        linintrup = ones(1,length(linx));
    end
    [f{i}, d{i}, r{i}] = filmthickness(CC,b{i},maxdistance,masksize,numorders,numaltmaps,linx,linperp,linintrup,s);
    F{i} = reshape(f{i},[s(1),s(2)]); % Min-distance film thickness
    D{i} = reshape(d{i},[s(1),s(2)]); % Distance
    R{i} = uint8(reshape(r{i},[s(1),s(2),3])); % Reconstructed image

    if sampcircrad>0
        A{i}(repmat(ringindex{i},1,1,3)) = 255;
        R{i}(repmat(outofcircindex{i},1,1,3)) = 0;
        F{i}(outofcircindex{i}) = NaN;
        D{i}(outofcircindex{i}) = 0;
        G{i} = F{i}(~outofcircindex{i});
        Favg(i) = mean(G{i}(G{i}~=0),'all');
    end
end

function [f,d,r] = filmthickness(C,RGB,maxdist,masksize,numorders,numaltmaps,linx,linperp,linintrup,sz) % these operations are not yet vectorised for speed
    diffkernel = [-2,1,1]; % This one works pretty well!
%     diffkernel = [1,-2,1];
    % assume the dimension-1 indices of C correspond to unit film thickness increments
    L2dist = sum((bsxfun(@minus,reshape(C(:,2:4),[size(C,1),1,3]),reshape(RGB(linx,:),[1,length(linx),3]))).^2,3); % sum of squares of colour diffs
    FF(1:numaltmaps) = {zeros(sz)};
    DD(1:numaltmaps) = {zeros(sz)};
    RR(1:numaltmaps) = {zeros([sz,3])};
    for j = 1:length(linx)
        [aa,bb] = ind2sub(sz,linx(j));
        for o = 1:numorders
            [dd(o),ii] = min(L2dist(:,j),[],1);
            Cimin = max(1,ii-masksize);
            Cimax = min(size(C,1),ii+masksize);
            L2dist(Cimin:Cimax,j) = Inf;
            iio(o) = ii;
        end
        if (j<6)
            for o = 1:numaltmaps
                FF{o}(linx(j)) = C(iio(o),1);
                DD{o}(linx(j)) = sqrt(dd(o));
                RR{o}(aa,bb,:) = C(iio(o),2:4);
            end
        else
            % Choose which o-value will minimise the gradient for each numorders
            % images we're creating
            for o = 1:numaltmaps
                dperp = zeros(1,numorders);
%                 dlin = zeros(1,numorders);
                fd = FF{o}(linperp(j,:));
                for oo = 1:numorders
                    fd(1) = C(iio(oo),1);
                    dperp(oo) = abs(sum(diffkernel.*fd));
                    if (~(linintrup(j)||linintrup(j-1)))
                        fd(2:3) = FF{o}(linx([j-1,j-2]));
                        dlin(oo) = abs(sum(diffkernel.*fd));
                    end
                end
                if (~(linintrup(j)||linintrup(j-1)))
                    [~,ooi] = min(dlin);
                else
                    [~,ooi] = min(dperp);
                end
                FF{o}(linx(j)) = C(iio(ooi),1);
                DD{o}(linx(j)) = sqrt(dd(ooi));
                RR{o}(aa,bb,:) = C(iio(ooi),2:4);
            end
        end
    end
    
%     Cimin = zeros(1,numorders);
%     Cimax = zeros(1,numorders);
%     % just do the first pixel
%     L2dist = sum((bsxfun(@minus,reshape(C(:,2:4),[size(C,1),1,3]),reshape(RGB(linx(1),:),[1,1,3]))).^2,3); % sum of squares of colour diffs
%     iii = true(1,size(C,1));
%     [dd,ii] = min(L2dist,[],1);
%     FF(1:numorders) = {zeros(sz)};
%     DD(1:numorders) = {zeros(sz)};
%     RR(1:numorders) = {zeros([sz,3])};
%     FF{1}(linx(1)) = C(ii,1);
%     DD{1}(linx(1)) = sqrt(dd);
%     [aa,bb] = ind2sub(sz,linx(1));
%     RR{1}(aa,bb,:) = C(ii,2:4);
%     Cimin(1) = max(1,ii-masksize);
%     Cimax(1) = min(size(C,1),ii+masksize);
%     for o = 2:numorders
%         iii(Cimin(o-1):Cimax(o-1)) = false;
%         L2dist(~iii) = Inf;
%         [dd,ii] = min(L2dist,[],1);
%         DD{o}(linx(1)) = sqrt(dd); % square root (L2 norm)
%         FF{o}(linx(1)) = C(ii,1);
%         RR{o}(aa,bb,:) = C(ii,2:4);
%         Cimin(o) = max(1,ii-masksize);
%         Cimax(o) = min(size(C,1),ii+masksize);
%     end
%     for j = 2:length(linx)
%         for o = 1:numorders
%             L2dist = sum((bsxfun(@minus,reshape(C(Cimin(o):Cimax(o),2:4),[Cimax(o)-Cimin(o)+1,1,3]),reshape(RGB(linx(j),:),[1,1,3]))).^2,3); % sum of squares of colour diffs
%             [dd,ii] = min(L2dist,[],1);
%             iio = ii+Cimin(o)-1;
%             DD{o}(linx(j)) = sqrt(dd); % square root (L2 norm)
%             FF{o}(linx(j)) = C(iio,1);
%             [aa,bb] = ind2sub(sz,linx(j));
%             RR{o}(aa,bb,:) = C(iio,2:4);
%             Cimin(o) = max(1,iio-masksize);
%             Cimax(o) = min(size(C,1),iio+masksize);
%         end
%         % Order swapping
%         dperp = zeros(1,numorders);
%         if (j>5) % after this point linperp only refers to pixels we have already processed
%             for o = 1:numorders
%                 fd = FF{o}(linperp(j,:));
%                 for oo = 1:numorders
%                     fd(1) = FF{oo}(linperp(j,1));
%                     dperp(oo) = abs(sum(diffkernel.*fd));
%                 end
%                 [~,ooi] = min(dperp);
%                 if ooi~=o
%                     swp = FF{ooi}(linx(j));
%                     FF{ooi}(linx(j)) = FF{o}(linx(j));
%                     FF{o}(linx(j)) = swp;
%                     swp = DD{ooi}(linx(j));
%                     DD{ooi}(linx(j)) = DD{o}(linx(j));
%                     DD{o}(linx(j)) = swp;
%                     swpv = RR{ooi}(aa,bb,:);
%                     RR{ooi}(aa,bb,:) = RR{o}(aa,bb,:);
%                     RR{o}(aa,bb,:) = swpv;
%                     swp = Cimin(ooi);
%                     Cimin(ooi) = Cimin(o);
%                     Cimin(o) = swp;
%                     swp = Cimax(ooi);
%                     Cimax(ooi) = Cimax(o);
%                     Cimax(o) = swp;
%                 end
%             end
%         end
%     end
    sd = cellfun(@(ce) sum(ce,'all'),DD);
    [~,oi] = min(sd);
    f = FF{oi};
    d = DD{oi};
    r = RR{oi};
    if maxdist>0
        f(d>maxdist) = 0;
        rgb(d>maxdist,:) = 0;
    end
end

end