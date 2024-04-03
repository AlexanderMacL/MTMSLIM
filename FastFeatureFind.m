% Function to find symmetries in an image by 1D convolution with optional
% averaging
% 
% I is an image of size m x n x 3 (RGB) or m x n (greyscale)
% avgblksz is a number specifying the number of elements of I that will be
%   averaged for each element in the autocorrelation, per dimension
% axis (optional) is a vector containing the axes of symmetry to be found,
%   if empty or not present, axes 1 and 2 are used

function [centers,radius] = FastFeatureFind(I,axis,ellipse,avgblksz,avgblksz2)
    plotfigs = false; % there is another one in imconvolveaxis below
    if nargin<5
        avgblksz2 = 4;
        if nargin<4
            avgblksz = 16;
            if nargin<3
                ellipse = false;
                if nargin<2
                    axis = [1,2];
                end
            elseif ellipse && length(axis)<2
                warning('When ellipse is true both axes are required, setting axis to [1 2]');
                axis = [1,2];
            end
        end
    end
    s = 1; % smoothing factor for radius finding
    centers = zeros(1,length(axis));
    radius = 0;
    A = imresize(I,1/avgblksz2);
    if (length(size(I))>2)
        AA = imgradient(A(:,:,1))+imgradient(A(:,:,2))+imgradient(A(:,:,3));
    else
        AA = imgradient(A);
    end
    AA = AA./max(AA,[],'all');
    AG = imresize(AA,avgblksz2/avgblksz);
    AG = AG./max(AG,[],'all');
    if plotfigs
        figure()
        imshow(I);
        figure()
        imshow(AG);
    end
    for a = axis
        centers(axis==a) = imconvolveaxis(AG,a);
        % do a second convolution which does the same thing but limiting the
        % range of Q to just the region around the centre from the 1st conv
        % transformed to a less heavily averaged space
        if (avgblksz>1)
            centers(axis==a) = (imconvolveaxis(AA,a,(centers(axis==a)-1)*avgblksz/avgblksz2,(centers(axis==a)+1)*avgblksz/avgblksz2)-1)*avgblksz2+1;
        end
        if centers(axis==a)<1; centers = zeros(1,length(axis)); warning("Q matrix minimum is outside region of image gradient"); return; end
    end
    if length(axis)>1 && all(axis==1 | axis==2) % circle
        [X1,X2] = size(AA);
        c1 = (centers(axis==1)-1)/avgblksz2+1;
        c2 = (centers(axis==2)-1)/avgblksz2+1;
        if (~ellipse)
            rext1 = min([c1,X1-c1+1,c2,X2-c2+1]);
            R1 = AA(c1+(0:rext1-1),c2).*AA(c1-(0:rext1-1),c2);
            R2 = AA(c1,c2+(0:rext1-1)).*AA(c1,c2-(0:rext1-1));
            [~,r] = max(smooth(R1'.*R2,s));
            radius = round(r*avgblksz2);
            if (plotfigs)
                figure();
                plot(smooth(R1'.*R2,s),'k');
            end
        else
            rext1 = min([c1,X1-c1+1]);
            rext2 = min([c2,X2-c2+1]);
            R1 = AA(c1+(0:rext1-1),c2).*AA(c1-(0:rext1-1),c2);
            [~,r1] = max(smooth(R1,s));
            R2 = AA(c1,c2+(0:rext2-1)).*AA(c1,c2-(0:rext2-1));
            [~,r2] = max(smooth(R2,s));
            radius = round([r1,r2].*avgblksz2);    
        end
    else % we assume axis only has one element
        if (axis==2); A = A'; end
        X1 = size(AA,1);
        c1 = (centers-1)/avgblksz2+1;
        rext1 = min(c1,X1-c1+1);
        R = sum(AA(c1+(0:rext1-1),:),2).*sum(AA(c1-(0:rext1-1),:),2);
        [~,r] = max(smooth(R,s));
        radius = (r-1)*avgblksz2+1;
    end
end

function [c] = imconvolveaxis(AG,axis,startpx,endpx)
    plotfigs = false;
    if (axis == 2); AG = AG'; end
    [X1,X2] = size(AG);
    if (nargin<3)
        startpx = 1;
        endpx = 2*X1-1;
    elseif startpx<1 || endpx<1
        c = 0;
        return;
    else
        startpx = round((startpx-1)*(2*X1-1)/X1) + 1;
        endpx = round((endpx-1)*(2*X1-1)/X1) + 1;
    end
    Q1 = zeros(endpx+1-startpx,1);
    AGX1e = zeros(3*X1-2,X2);
    AGX1e(1:X1-1,:) = repmat(AG(1,:),X1-1,1);
    AGX1e(X1:2*X1-1,:) = AG;
    AGX1e(2*X1:3*X1-2,:) = repmat(AG(X1,:),X1-1,1);
    AGflipX1 = AG(X1:-1:1,:);
    for i=startpx:endpx
        Q1(i-startpx+1) = sum((AGflipX1-AGX1e(i:i+X1-1,:)).^2,'all');
    end
    if plotfigs
        figure()
        if (axis==2); imshow(AGX1e'); else; imshow(AGX1e); end
        figure()
        if (axis==2); imshow(repmat(Q1'./max(Q1,[],'all'),[20,1])); else; imshow(repmat(Q1./max(Q1,[],'all'),[1,20])); end
    end
    [~,ii1] = min(Q1);
%     c = ii1 - (X1 - 1);
    c = round((ii1-1+startpx)./(2*X1-2).*X1);
    if (c<1 || c>X1)
%         warning(['Q matrix minimum is outside region of image gradient in axis ' num2str(axis)]);
        c = 0;
        return;
    end
end