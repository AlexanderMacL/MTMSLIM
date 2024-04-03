% MTM SLIM Calculation

% Author: Alexander MacLaren
% Last updated: 01/08/2023

% Arguments:
%   testpath is a string or character vector containing the directory in
%       which the mtmd file is located, including trailing /
%   testname is a string or character vector containing the name of the
%       mtmd file without file extension
% All other arguments match those in SLIMwrap.m and are detailed there

% Return values:
%   A is a cell array containing RGB images showing contact
%   AA is a cell array containing cropped RGB images
%   F is a cell array containing contact film thickness maps
%   D is a cell array containing distance maps
%   R is a cell array containing reconstructed images
%   fnames is a cell array of strings containing the image filenames
%   Favg is a vector containing the average film thickness
%   mapi is a logical vector indicating the mapper steps
%   J is the datastructure extracted from the mtmd/etmd file

% For usage example see MTMSLIMexample.m

function [A,AA,F,D,R,fnames,Favg,mapi,J] = MTMSLIM(testpath, testname, ...
    maxdistance, minfilm, maxfilm, circlefitting, avgblocksize, ...
    sampcircrad, sampcirccent, calcurverefit, calcurveoffsetprescale, ...
    numorders, numaltmaps, masksize, tangentcorrection, plotfigures)

fnames = {};
Favg = [];
[A,AA,F,D,R] = deal(0);
B = [];

% extract filenames from mtmd file
J = mtmd2json([testpath, testname]);
mapi = cellfun(@(x) strcmp(x.stepType,'Mapper'),J.Steps);
blanki = cellfun(@(x) isempty(x.imgName), J.Steps(mapi));
if any(blanki) % remove mapper images never taken
    c=0;
    for i=1:length(mapi)
        if mapi(i)
            c=c+1;
            if blanki(c)
                mapi(i) = 0;
            end
        end
    end
end
c = sum(mapi);
if c==0; return; end

testname = testname(1:end-5);

imagepath = [testpath testname ' mapper images\'];
imagenames = cellfun(@(x) x.imgName,J.Steps(mapi),'UniformOutput',false);
calpath = [imagepath,testname,'-3D_SpacerCalibration.txt'];

[A,AA,F,D,R,Favg] = SLIMwrap(imagepath,imagenames,calpath,maxdistance,minfilm,maxfilm,circlefitting,avgblocksize,sampcircrad,sampcirccent,calcurverefit,calcurveoffsetprescale,numorders,numaltmaps,masksize,tangentcorrection,plotfigures);

end
