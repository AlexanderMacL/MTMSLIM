%% Function to convert *.mtmd or *.etmd file to JSON encoded data
%
% author: Alexander MacLaren
% revised: 22/11/2021
%
% Usage:
%   K = mtmd2json() - a file open dialog is provided to open *.mtmd file
%       and to save *.json file. To be used at the command line or as a
%       standalone script.
%   K = mtmd2json(infileloc) - the string or character vector infileloc
%       specifies the location of the input (mtmd) file - no output file is
%       written.
%   K = mtmd2json(infileloc, outfileloc) - the strings or character vectors
%       infileloc and outfileloc respectively contain the locations of the
%       desired input (mtmd) and output (JSON) files.
% 
% Return value:
%   K is a MATLAB data structure containing the contents of the mtmd file
%   as received by jsonencode()
%
% Notes:
%   ETM profile files are also supported - the filetypes are encoded
%       differently and the encoding is inferred from the input filepath
%       extension provided
%   The MTM file is expected to begin 0x02 0x00 0x00 0x00, and the expected
%       ETM file header is 0x01 0x00 0x00 0x00. If 'Warning: Unexpected
%       pattern in header' is thrown quoting one of these patterns,
%       consider changing the file extension as appropriate
%

function [K] = mtmd2json(varargin)

if (nargin==0)
    % open input file
    [flnm,pth,~] = uigetfile({'*.mtmd','MTM Data Files (*.mtmd)';'*.etmd','ETM Data Files (*.etmd)';'*.*','All Files (*.*)'});
    f = fopen([pth,flnm]);

elseif (nargin>2)
    error("Too many arguments (%d given)",nargin);

else 
    flnm = varargin{1};
    f = fopen(flnm);
end

mtm_etm = flnm(length(flnm)-3);
raw = fread(f, inf, 'uint8=>uint8');
data = raw';
fclose(f);

disp("Parsing file "+flnm+" ...");

c = uint64(length(data));
k = uint64(1);
i = 0;

if (c>=34)
    if (mtm_etm == 'm') % MTM
        dmtm = [2 0 0 0]==data(k:k+4-1);
        if (~all(dmtm))
            warning("Unexpected pattern " + num2str(data(k:k+4-1)) + " in header - is this an MTM data file?");
        end
    elseif (mtm_etm == 'e') % ETM
        detm = [1 0 0 0]==data(k:k+4-1);
        if (~all(detm))
            warning("Unexpected pattern " + num2str(data(k:k+4-1)) + " in header - is this an ETM data file?");
        end
    else
        error("File format of %s unsupported, supported formats are *.mtmd and *.etmd",flnm);
    end
    k = k + 4;
    [K.dataFilePath, ki] = parsestr(data, k); k = k + ki;
    [K.profileFilePath, ki] = parsestr(data, k); k = k + ki;
    [K.description, ki] = parsestr(data, k); k = k + ki;
    [K.lubeName, ki] = parsestr(data, k); k = k + ki;
    [K.comments, ki] = parsestr(data, k); k = k + ki;
    K.numSteps = typecast(data(k:k+4-1),'uint32'); k = k + 4;
    K.testStartTime = char(string(datetime(typecast(data(k:k+8-1), 'uint64'), 'ConvertFrom', '.net', 'TimeZone', 'Europe/London'),'dd/MM/yyyy HH:mm:ss','en_GB')); k = k + 8;
    [K.status, ki] = parsestr(data, k); k = k + ki;
    K.testEndTime = char(string(datetime(typecast(data(k:k+8-1), 'uint64'), 'ConvertFrom', '.net', 'TimeZone', 'Europe/London'),'dd/MM/yyyy HH:mm:ss','en_GB')); k = k + 8;
    K.numStepsCompleted = typecast(data(k:k+4-1),'uint32'); k = k + 4;
else
    error("File shorter than expected header");
end

while (k<c)
    stephead = typecast(data(k:k+8-1),'uint64');
    k = k + 8;
    i = i + 1;
    switch (stephead)
        case {hex2dec('0000000200000000'), hex2dec('0000000300000000')} % Traction
            disp ("Traction step "+num2str(i)+" of "+num2str(K.numStepsCompleted));
            K.Steps{i}.stepType = 'Traction';
        case {hex2dec('0000000200000001'), hex2dec('0000000300000001')} % Stribeck
            disp ("Stribeck step "+num2str(i)+" of "+num2str(K.numStepsCompleted));
            K.Steps{i}.stepType = 'Stribeck';
        case {hex2dec('0000000200000002'), hex2dec('0000000300000002')} % Timed
            disp ("Timed step "+num2str(i)+" of "+num2str(K.numStepsCompleted));
            K.Steps{i}.stepType = 'Timed';
        case {hex2dec('0000000200000008'), hex2dec('0000000300000003')} % Mapper
            disp ("Mapper step "+num2str(i)+" of "+num2str(K.numStepsCompleted));
            K.Steps{i}.stepType = 'Mapper';
        otherwise
            error("Step header "+num2str(typecast(data(k-8:k-1),'uint32'))+" at byte "+num2str(k-8)+" unrecognised")
    end
    [K.Steps{i}.stepName, ki] = parsestr(data, k); k = k + ki;
    K.Steps{i}.stepStartTime = char(string(datetime(typecast(data(k:k+8-1), 'uint64'), 'ConvertFrom', '.net', 'TimeZone', 'Europe/London'),'dd/MM/yyyy HH:mm:ss','en_GB')); k = k + 8;
    if(strcmp(K.Steps{i}.stepType, 'Timed'))
        K.Steps{i}.stepDurationSeconds = round(double(typecast(data(k:k+8-1), 'uint64'))/1e7,3);
    end
    k = k + 8;
    K.Steps{i}.tractionForceZero = round(typecast(data(k:k+8-1), 'double'),3); k = k + 8;
    K.Steps{i}.tractionForceZeroMeasured = data(k)==1; k = k + 1;
    if (strcmp(K.Steps{i}.stepType, 'Mapper'))
        if (mtm_etm == 'e') % ETM
            k = k + 15;
        else
            k = k + 13; % skip the data that isn't relevant to Mapper
        end
        if (typecast(data(k:k+4-1),'uint32')~=2)
            warning("Unexpected Mapper header sequence %d at byte %d", (typecast(data(k:k+4-1),'uint32')), k);
        end
        k = k + 4;
        K.Steps{i}.windowLoad = round(typecast(data(k:k+8-1),'double'),3); k = k + 8;
        [K.Steps{i}.imgName, ki] = parsestr(data, k); k = k + ki;
        k = k + 1;
        if(data(k-1)==1) % If first mapper image (ZERO)
            K.Steps{i}.selectorCentre(1) = typecast(data(k:k+4-1),'uint32'); k = k + 4;
            K.Steps{i}.selectorCentre(2) = typecast(data(k:k+4-1),'uint32'); k = k + 4;
            K.Steps{i}.selectorRadius = typecast(data(k:k+4-1),'uint32'); k = k + 4;
            K.Steps{i}.spacerLayerThickness = round(typecast(data(k:k+8-1),'double'),3); k = k + 8;
        end
    else
        K.Steps{i}.trackRadiusCalibrationEn = data(k)==1; k = k + 1;
        K.Steps{i}.trackRadiusCalibrationSuccessful = data(k)==1; k = k + 1;
        k = k + 3; % Not sure what these are for
        K.Steps{i}.discTrackRadius = round(typecast(data(k:k+8-1),'double'),3); k = k + 8;
        if (mtm_etm == 'e') % ETM
            K.Steps{i}.trackRadiusCompensationEn = data(k)==1; k = k + 1;
            K.Steps{i}.skewCompensationEn = data(k)==1; k = k + 1;
        end
        if (typecast(data(k:k+4-1),'uint32')~=1)
            warning("Unexpected datapoints header sequence %d at byte %d", (typecast(data(k:k+4-1),'uint32')), k);
        end
        k = k + 4;
        K.Steps{i}.numDatapoints = typecast(data(k:k+4-1),'uint32'); k = k + 4;
        for j = 1:K.Steps{i}.numDatapoints
            % K.Steps{i}.dataTimestamp(j) = string(datetime(typecast(data(k:k+8-1), 'uint64'), 'ConvertFrom', '.net', 'TimeZone', 'Europe/London'),'dd/MM/yyyy HH:mm:ss','en_GB');
            k = k + 8;
            K.Steps{i}.secondsElapsed(j) = round(typecast(data(k:k+8-1),'double'),3); k = k + 8;
            K.Steps{i}.lubeTemperature(j) = round(typecast(data(k:k+8-1),'double'),2); k = k + 8;
            K.Steps{i}.potTemperature(j) = round(typecast(data(k:k+8-1),'double'),2); k = k + 8;
            K.Steps{i}.ballLoad(j) = round(typecast(data(k:k+8-1),'double'),3); k = k + 8;
            K.Steps{i}.wear(j) = round(typecast(data(k:k+8-1),'double'),1); k = k + 8;
            if (mtm_etm ~= 'e') % ETM has no ECR
                % I think this is ECR but not totally sure
                % ignoring these for now
                % K.Steps{i}.ECR(j) = typecast(data(k:k+8-1),'double');
                k = k + 8;
            end
            if (strcmp(K.Steps{i}.stepType, 'Traction') || strcmp(K.Steps{i}.stepType, 'Stribeck'))
                K.Steps{i}.ballSpeed1(j) = round(typecast(data(k:k+8-1),'double'),1); k = k + 8;
                K.Steps{i}.ballSpeed2(j) = round(typecast(data(k:k+8-1),'double'),1); k = k + 8;
                K.Steps{i}.discSpeed1(j) = round(typecast(data(k:k+8-1),'double'),1); k = k + 8;
                K.Steps{i}.discSpeed2(j) = round(typecast(data(k:k+8-1),'double'),1); k = k + 8;
                K.Steps{i}.rollingSpeed(j) = round(typecast(data(k:k+8-1),'double'),2); k = k + 8;
                K.Steps{i}.slidingSpeed(j) = round(typecast(data(k:k+8-1),'double'),2); k = k + 8;
                K.Steps{i}.SRR(j) = round(typecast(data(k:k+8-1),'double'),2); k = k + 8;
                K.Steps{i}.TF1(j) = round(typecast(data(k:k+8-1),'double'),2); k = k + 8;
                K.Steps{i}.TF2(j) = round(typecast(data(k:k+8-1),'double'),2); k = k + 8;
            elseif (strcmp(K.Steps{i}.stepType, 'Timed'))
                K.Steps{i}.ballSpeed(j) = round(typecast(data(k:k+8-1),'double'),1); k = k + 8;
                K.Steps{i}.discSpeed(j) = round(typecast(data(k:k+8-1),'double'),1); k = k + 8;
                K.Steps{i}.rollingSpeed(j) = round(typecast(data(k:k+8-1),'double'),2); k = k + 8;
                K.Steps{i}.slidingSpeed(j) = round(typecast(data(k:k+8-1),'double'),2); k = k + 8;
                K.Steps{i}.SRR(j) = round(typecast(data(k:k+8-1),'double'),2); k = k + 8;
            end
            K.Steps{i}.tractionForce(j) = round(typecast(data(k:k+8-1),'double'),3); k = k + 8;
            K.Steps{i}.tractionCoefficient(j) = round(typecast(data(k:k+8-1),'double'),4); k = k + 8;
        end
    end
end

if (nargin==0 || nargin==2)
    str = jsonencode(K); %, 'PrettyPrint', true);
    % PrettyPrint option not working until R2021a so this is the poor man's JSON beautifier (no indentation)
    str = strrep(str, ',"', sprintf(',\n"'));
    str = strrep(str, '":', '": ');
    str = strrep(str, '{', sprintf('\n{\n'));
    str = strrep(str, '}', sprintf('\n}\n'));
    str = strrep(str, sprintf('}\n,'), '},');

    if (nargin==2)
        flnm = varargin{2};
        f = fopen(flnm,'w');
    elseif (nargin==0)
        [flnm,pth,~] = uiputfile({'*.json','JavaScript Object Notation Files (*.json)';'*.*','All Files (*.*)'},'Save File',flnm(1:length(flnm)-5));
        f = fopen([pth,flnm],'w');
    end

    % write JSON string to file
    fwrite(f, str(2:end));
    fclose(f);

end
end

function [str, ki] = parsestr(data, k)
    ki = 1;
    str = '';
    c = uint64(data(k));
    if (c==0)
        return;
    else
        ci = 1;
        while (bitand(uint64(data(k)),128)~=0)
            c = c + bitshift(1,7*ci)*(uint64(data(k+1))-1); % if more than 127 chars in descriptor, account for extra byte
            ki = ki + 1;
            k = k + 1;
            ci = ci + 1;
        end
    end
    str = char(data(k+1:k+c));
    ki = c + ki;
end
