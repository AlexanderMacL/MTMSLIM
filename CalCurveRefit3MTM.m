function [CC] = CalCurveRefit3MTM(C,CC,betaR,betaG,betaB,b,minfilm,maxfilm,calcurveoffsetprescale,plotfigures,numpixels)
    plot3Dgraph(C,CC,[],'Initial Fit',plotfigures)
    disp('Preparing to refit calibration curve');
    nonzeroindices = ~(b(:,1)==0 & b(:,2)==0 & b(:,3)==0);
    nonzeroindices = nonzeroindices & (isfinite(b(:,1)) & isfinite(b(:,2)) & isfinite(b(:,3)));
    RGB = double(b);
    RGBnz = RGB(nonzeroindices,:);
    if (numpixels>0)
        % reduce number of pixels in image
        skip = round(size(RGBnz,1)/numpixels);
        if (skip<1); skip=1; end
        RGBnz = RGBnz(1:skip:end,:);
    end
    
    [f] = filmthickness(CC,RGBnz);
    plot3Dgraph(C,CC,RGBnz,'Before Refit',plotfigures)
    plot2Dgraph(C,CC,f,RGBnz,'Before Refit',plotfigures)

    mx = max(RGBnz,[],1);
    mn = min(RGBnz,[],1);
    mxc = max(CC(:,2:4),[],1);
    mnc = min(CC(:,2:4),[],1);
    % newcolour = oldcolour*(C1+C3h) + C2
%     beta([1,2,3]) = (log((mx-mn)./(mxc-mnc))+[betaR(3),betaG(3),betaB(3)]); % gradient
    beta([1,2,3]) = log((mx-mn)./(mxc-mnc))+[betaR(3),betaG(3),betaB(3)];
%     beta([4,5,6]) = (mnc.*mx - mn.*mxc)./(mnc-mxc); % offset
    beta([4,5,6]) = mn - (mnc+calcurveoffsetprescale).*(mx - mn)./(mxc-mnc); % offset
%     beta([4,5,6]) = [0,0,0];
%     beta = [1,1,1,0,0,0];
    beta([7,8,9]) = [betaR(4),betaG(4),betaB(4)];
    
    numiter = 3;
%     CC = C;
    % do initial cal curve shift
    CC(:,2:4) = film2colour(CC(:,1),betaR,betaG,betaB,beta);
    
    plot3Dgraph([],CC,RGBnz,'Initial Refit',plotfigures);
    
    for i = 1:numiter
        disp(['Refitting calibration curve iteration ' num2str(i) ' of ' num2str(numiter)]);
        % calculate film thicknesses for this calibration curve
        [f] = filmthickness(CC,RGBnz);
        % optimise fitting parameters maintaining prior film thickness calc
        opts = statset('nlinfit');
        opts.RobustWgtFun = 'huber';
        beta = nlinfit(f,zeros(size(RGBnz,1),1),@(b, h) sum((film2colour(h,betaR,betaG,betaB,b)-RGBnz).^2,2), beta);
        % recalculate cal curve
        CC(:,2:4) = film2colour(CC(:,1),betaR,betaG,betaB,beta);
    end
    % crop cal curve
    if maxfilm ~= 0; CC = CC(CC(:,1)<=maxfilm,:); end
    if minfilm ~= 0; CC = CC(CC(:,1)>=minfilm,:); end
    
%     % do final film thickness calc
    [f] = filmthickness(CC,RGBnz);
%     [f] = filmthickness(CC,RGBnz);
     
    plot3Dgraph([],CC,RGBnz,'After Refit',plotfigures)
    plot2Dgraph(C,CC,f,RGBnz,'After Refit',plotfigures);
end

function [cc] = film2colour(h,betaR,betaG,betaB,beta) % just a placeholder for however you want to lookup the calibration
    betaR(3:4) = beta([1,7]);
    betaG(3:4) = beta([2,8]);
    betaB(3:4) = beta([3,9]);
    r = expexpdecsin(betaR,h);
    g = expexpdecsin(betaG,h);
    b = expexpdecsin(betaB,h);
    cc = [r,g,b]+beta(4:6);
%     cc = [r,g,b].*(beta(1:3) + h.*beta(7:9))+beta(4:6); % gradient h-dependent
%     cc = [r,g,b].*(beta(1:3))+beta(4:6)+h.*beta(7:9); % offset h-dependent
end

function [f] = filmthickness(C,RGB) % these operations are vectorised for speed
    L2dist = sum((bsxfun(@minus,reshape(C(:,2:4),[size(C,1),1,3]),reshape(RGB,[1,size(RGB,1),3]))).^2,3); % sum of squares of colour diffs
    [~,ii] = min(L2dist,[],1);
    f = C(ii,1);
end

function plot3Dgraph(C,CC,RGBnz,ttl,plotfigures)
    if (plotfigures)
        figure()
        if (~isempty(C))
            idx = 1:3:size(C,1);
            scatter3(C(idx,2),C(idx,3),C(idx,4),ones(1,length(idx))*30,C(idx,2:4)./255,'s','filled','DisplayName','Cal');
            hold on
        end
        scatter3(CC(:,2),CC(:,3),CC(:,4),ones(1,size(CC,1))*20,CC(:,2:4)./255,'o','filled','DisplayName','Cal fit');
        hold on
        if (~isempty(RGBnz))
            scatter3(RGBnz(:,1),RGBnz(:,2),RGBnz(:,3),ones(1,size(RGBnz,1))*2,RGBnz./255,'filled','DisplayName','pixels');
        end
        set(gca,'FontSize',16);
        xlabel('Red');
        ylabel('Green');
        zlabel('Blue');
        axis equal
        view(135,30); %40,30 %-20,18
        grid on
        xticks([0:50:250]);
        yticks([0:50:250]);
        zticks([0:50:250]);
        xlim([0,260]);
        ylim([0,260]);
        zlim([0,260]);
        ll=legend('Location','ne');
        ll.BoxFace.ColorType='truecoloralpha';
        ll.BoxFace.ColorData=uint8(255*[1 1 1 0.7]');
        title(ttl)
        drawnow;
    end
end

function plot2Dgraph(C,CC,f,RGB,ttl,plotfigures)
    if plotfigures
        figure()
        colororder([1,0,0;0,1,0;0,0,1]);
        plot(C(:,1),C(:,2:4),'s','MarkerSize',4,'MarkerFaceColor','auto');
        hold on
        plot(CC(:,1),CC(:,2:4));
        scatter(f,RGB,4,'.');
        set(gca,'FontSize',16);
        leg = {'cal R','cal G','cal B','fit R','fit G','fit B','image R','image G','image B'};
        xlabel('Film thickness (nm)')
        ylabel('Intensity (LSB)')
        ylim([0,255]);
        ll=legend(leg,'NumColumns',length(leg)/3,'Location','ne');
        ll.BoxFace.ColorType='truecoloralpha';
        ll.BoxFace.ColorData=uint8(255*[1 1 1 0.7]');
        title(ttl);
    end
end

