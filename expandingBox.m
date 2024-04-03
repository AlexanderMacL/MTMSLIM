% function to calculate linear indices in order for expanding box

% 'Expanding box' - to locate pixels on box edge given centre pixel (a,b)
% consider the ith iteration, where 0 is where there is only 1 pixel
% 0 start at (a,b-i), 1 move i to (a+i,b-i)
% 2 move 2i to (a+i,b+i), 3 move 2i to (a-i,b+i),
% 4 move 2i to (a-i,b-i), 5 move i-1 to (a-i+1,b-i)
%   b ->
% a 4 4 4 4 3
% | 5 4 4 3 3 
% v 0 0 . 3 3 
%   1 1 2 2 3
%   1 2 2 2 2
%
% given co-ordinates (m,n) from above, use for gradient:
% 0 (m,n), (m,n+1) and (m,n+2)
% 1 etc: (m,n) and the last two from above or the last two previously

% actually it would be better if the traverse changed rotational direction and incremented i every time
% we reached a 'boundary' (either edge of image or the end of this
% rotation)

% 'direction' just means that 1, 3 and 5 go the other way
% and we shall assume the image is square and not allow the traverse to
% bifurcate?

% linperp will have a 2nd dimension of length 3 - these are the elements to
% use to take a perpendicular 2nd derivative with direction 'into' the
% square

% THE WHOLE LINPERP THING ASSUMES WE START INSIDE THE CIRCLE AND CURRENTLY
% DOESN'T HANDLE CORNERS PROPERLY

function [linx,linperp,linintrup] = expandingBox(shape,startpx,outofcircindex)
    linx = zeros(1,shape(1)*shape(2));
    linperp = zeros(shape(1)*shape(2),3);
    linintrup = false(1,shape(1)*shape(2));
    if (length(shape)>2); warning(['Extra dimensions ignored for shape vector >2D (' num2str(shape) ')']); end
    if (length(shape)>2); warning(['Extra dimensions ignored for startpx >2D (' num2str(startpx) ')']); end
    a = startpx(1);
    b = startpx(2);
    linx(1) = sub2ind(shape,a,b);
    ii = 2;
    dir = 1;
    for i=1:max([startpx, shape-startpx])
        for j = (dir-1)/2:dir:i*dir
            if a+j>shape(1) || a+j<1 || b-i<1; linintrup(ii) = true; continue; end
            linx(ii) = sub2ind(shape,a+j,b-i);
            if (outofcircindex(linx(ii-1))); linintrup(ii) = true; end
            if (j==i*dir)
                linperp(ii,1:3) = sub2ind(shape,a+j-[0,1,2]*dir,b-i+[0,1,2]);
            else
                linperp(ii,1:3) = sub2ind(shape,a+j+[0,0,0],b-i+[0,1,2]);
            end
            ii = ii + 1;
        end
        for j = 1-i:i
            if a+i*dir>shape(1) || a+i*dir<1 || b+j<1 || b+j>shape(2); linintrup(ii) = true; continue; end
            linx(ii) = sub2ind(shape,a+i*dir,b+j);
            if (outofcircindex(linx(ii-1))); linintrup(ii) = true; end
            if (j==i)
                linperp(ii,1:3) = sub2ind(shape,(a+(i-[0,1,2])*dir),b+j-[0,1,2]);
            else
                linperp(ii,1:3) = sub2ind(shape,(a+(i-[0,1,2])*dir),b+j+[0,0,0]);
            end
            ii = ii + 1;
        end
        for j = (i-1)*dir:-1*dir:-i*dir
            if b+i>shape(2) || a+j<1 || a+j>shape(1); linintrup(ii) = true; continue; end
            linx(ii) = sub2ind(shape,a+j,b+i);
            if (outofcircindex(linx(ii-1))); linintrup(ii) = true; end
            if (j==-i*dir)
                linperp(ii,1:3) = sub2ind(shape,a+j+[0,1,2]*dir,b+i-[0,1,2]);
            else
                linperp(ii,1:3) = sub2ind(shape,a+j+[0,0,0],b+i-[0,1,2]);
            end
            ii = ii + 1;
        end
        for j = i-1:-1:-i
            if a-i*dir>shape(1) || a-i*dir<1 || b+j<1 || b+j>shape(2); linintrup(ii) = true; continue; end
            linx(ii) = sub2ind(shape,a-i*dir,b+j);
            if (outofcircindex(linx(ii-1))); linintrup(ii) = true; end
            if (j==-i)
                linperp(ii,1:3) = sub2ind(shape,a-(i-[0,1,2])*dir,b+j+[0,1,2]);
            else
                linperp(ii,1:3) = sub2ind(shape,a-(i-[0,1,2])*dir,b+j+[0,0,0]);
            end
            ii = ii + 1;
        end
        for j = (1-i)*dir:dir:(dir+1)/-2
            if b-i<1 || a+j<1 || a+j>shape(1); linintrup(ii) = true; continue; end
            linx(ii) = sub2ind(shape,a+j,b-i);
            if (outofcircindex(linx(ii-1))); linintrup(ii) = true; end
            linperp(ii,1:3) = sub2ind(shape,a+j+[0,0,0],b-i+[0,1,2]);
            ii = ii + 1;
        end
        dir = dir*-1;
%         for j = 0:i
%             if a+j>shape(1) || b-i<1; continue; end
%             linx(ii) = sub2ind(shape,a+j,b-i);
%             ii = ii + 1;
%         end
%         for j = 1-i:i
%             if a+i>shape(1) || b+j<1 || b+j>shape(2); continue; end
%             linx(ii) = sub2ind(shape,a+i,b+j);
%             ii = ii + 1;
%         end
%         for j = i-1:-1:-i
%             if b+i>shape(2) || a+j<1 || a+j>shape(1); continue; end
%             linx(ii) = sub2ind(shape,a+j,b+i);
%             ii = ii + 1;
%         end
%         for j = i-1:-1:-i
%             if a-i<1 || b+j<1 || b+j>shape(2); continue; end
%             linx(ii) = sub2ind(shape,a-i,b+j);
%             ii = ii + 1;
%         end
%         for j = 1-i:-1
%             if b-i<1 || a+j<1; continue; end
%             linx(ii) = sub2ind(shape,a+j,b-i);
%             ii = ii + 1;
%         end
    end
    if ii<shape(1)*shape(2); warning('missing indices'); end
    % remove bits out of circle
    linintrup = linintrup(~outofcircindex(linx));
    linperp = linperp(~outofcircindex(linx),:);
    linx = linx(~outofcircindex(linx));
end