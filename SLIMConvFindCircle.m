% Function to convolutionally find SLIM circle

function [circrad, circloc] = SLIMConvFindCircle(I, avgblksz)

    if nargin<2
        avgblksz = 16; % average block size
    end
    
    s = 5; % smoothing for radius finding
    minr = 64; % minimum radius in pixels (original image)
    avgr = 4;
    % ratio of 1st convolution block size to 2nd convolution block size
    % (higher is more computationally expensive but improves accuracy of
    % both location and radius detection)
    avbs = avgblksz/avgr;
    
    circrad = 0;
    circloc = [0,0];
    
%     AAa = imgaverage(I,avbs);
    AAa = imresize(I,1/avbs);
%     A = imgaverage(AAa,avgr);
    A = imresize(AAa,1/avgr);
    AG = imgradient(A(:,:,1))+imgradient(A(:,:,2))+imgradient(A(:,:,3));
    AG = AG./max(AG,[],'all');

    [Y,X,~] = size(AG);
    Q = imgconvolve90(AG);
    [~,xlin] = min(Q,[],'all','linear');
    [y,x,~] = ind2sub(size(Q),xlin);
    if (y < X/2-1) || (y > X/2-1+Y-1) || (x < Y/2-1) || (x > Y/2-1+X-1)
        warning(['Q matrix minimum (' num2str(y) ',' num2str(x), ') is outside region of original image']);
        return
    end

%     a = x-1-(X+Y-2)/2;
%     b = y-1-(X+Y-2)/2;
%     u = floor((a-b)/2+X/2+1);
%     v = floor(Y/2+1+(a+b)/2);

    % Second convolution (correction)
    AAG = imgradient(AAa(:,:,1))+imgradient(AAa(:,:,2))+imgradient(AAa(:,:,3));
    AAG = AAG./max(AAG,[],'all');
    [YY,XX,~] = size(AAG);
    QQ = imgconvolve90(AAG,(x-1)*avgr,(x+1)*avgr,(y-1)*avgr,(y+1)*avgr);
    [~,xxlin] = min(QQ,[],'all','linear');
    [yy,xx,~] = ind2sub(size(QQ),xxlin);

    aa = xx-1-(XX+YY-2)/2;
    bb = yy-1-(XX+YY-2)/2;
    uu = floor((aa-bb)/2+XX/2+1);
    vv = floor(YY/2+1+(aa+bb)/2);
    
    if (uu < 1 || vv < 1 || uu > XX || vv > yy)
        warning(['Second convolution moved centre (' num2str(uu) ',' num2str(vv) ') outside image']);
        return
    end

    Aa = uu;
    Bb = vv;
    ofs = floor(minr/avbs);
    Vr = smooth(AAG(Bb,Aa:XX),s);
    Vl = smooth(flip(AAG(Bb,1:Aa)),s);
    Vt = smooth(flip(AAG(1:Bb,Aa)),s);
    Vb = smooth(AAG(Bb:YY,Aa),s);
    ml = max([length(Vr), length(Vl), length(Vt), length(Vb)]);
    Vr = padarray(Vr,ml-length(Vr),'post');
    Vl = padarray(Vl,ml-length(Vl),'post');
    Vt = padarray(Vt,ml-length(Vt),'post');
    Vb = padarray(Vb,ml-length(Vb),'post');
    Vff = Vr.*Vl.*Vt.*Vb;
    Vff(1:ofs) = 0;

    [~,iff] = max(Vff);
    circrad = (iff-1)*avbs;
    circloc = [Aa,Bb].*avbs;
end


function [Iavg] = imgaverage(I,n)
    crosshairs = floor(size(I,[1,2])./2); % zero-indexed co-ordinates
    avgex = floor(crosshairs(2)/n);
    avgey = floor(crosshairs(1)/n);
    Iavg = zeros(2*avgey,2*avgex,size(I,3));
    for p = -avgey:avgey-1
        for q = -avgex:avgex-1
            Iavg(p+avgey+1,q+avgex+1,:) = mean(I(crosshairs(1)+p*n+1:crosshairs(1)+(p+1)*n, crosshairs(2)+q*n+1:crosshairs(2)+(q+1)*n, :),[1,2]);
        end
    end
    Iavg = uint8(Iavg);
end

function [Q] = imgconvolve90(I,startx,finishx,starty,finishy)
    [Y,X,~] = size(I);
    
    if nargin<5
        finishx = X+Y-2;
        startx = 1;
        finishy = X+Y-2;
        starty = 1;
    end
    
    % Extend image
    E(X:X+Y-1,Y:Y+X-1) = I;

    E(1:X-1,Y:Y+X-1) = repmat(I(1,:),X-1,1); % top edge
    E(X+Y:2*X-2+Y,Y:Y+X-1) = repmat(I(Y,:),X-1,1); % bottom edge
    E(X:X+Y-1,1:Y-1) = repmat(I(:,1),1,Y-1); % left edge
    E(X:X+Y-1,Y+X:X+2*Y-2) = repmat(I(:,X),1,Y-1); % right edge
    E(1:X-1,1:Y-1) = repmat(I(1,1),X-1,Y-1); % top-left corner
    E(1:X-1,Y+X:X+2*Y-2) = repmat(I(1,X),X-1,Y-1); % top-right corner
    E(X+Y:2*X-2+Y,1:Y-1) = repmat(I(Y,1),X-1,Y-1); % bottom-left corner
    E(X+Y:2*X-2+Y,Y+X:X+2*Y-2) = repmat(I(Y,X),X-1,Y-1); % bottom-right corner
    B = imrotate(I,-90);
    Q = ones(Y+X-2,X+Y-2);
    M = true(Y+X-2,X+Y-2);
    M(starty:finishy,startx:finishx) = false;
    for p = starty:finishy
        for q = startx:finishx
            Q(p,q) = sum((abs(E(p:p+X-1,q:q+Y-1)-B(:,:))).^2,[1,2]);
        end
    end
    Q = Q./max(Q,[],'all');
    Q(M) = NaN;
end