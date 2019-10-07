function values = get_goalc_parameters( filePath, objectStr, variableArray )
%GET_GOALC_PARAMETERS  extract specified parameters values from goalc.txt file
%
%   GET_GOALC_PARAMETERS( filePath, objectStr, variableArray )
% 
%   Output:
%       values              - values corresponding to specified variables -
%                             can be single or multiple floats, or string
%
%   Input:
%       filePath            - file path to input goalc.txt file
%       objectStr           - string specifying goalc object to parse
%       variableArray       - array of strings specifying variables within goalc object
%
%   Example:
%       values = mrecon_pcmr_parsegoalc( <../goalc.txt>, 'GR `ex', {'slope', 'dur'} )
%
%   Usage:
%       Use with Josh's mrecon function: mrecon_writegoalc2txt.m

% TAR (t.roberts@kcl.ac.uk)


%% TO DO:
%- error if object or variables not found
%- some goalc objects have multiple composites - currently only accesses
%  variables within first composite. Might need flexibility to access
%  further than 1 composites deep.
%- speed-up: loop within an object block rather than looping through the goalc file every time


%% Find index of line containing object string

fid = fopen(filePath);
tline = fgetl(fid);
tlineCounter = 1;

% loop line-by-line until find objectStr
while ischar(tline)
    if strfind(tline, objectStr)
        break;
    end
    % Read next line in fid
    tline = fgetl(fid);
    tlineCounter = tlineCounter + 1;
end
fclose(fid);

% return line Idx of objectStr
objectLine = tlineCounter;


%% check if variableArray is singular
if ischar(variableArray)
    variableArray = {variableArray};
end


%% Find indices of lines containing variable strings

for ii = 1:numel(variableArray)
    
    clear varCurrent
    varCurrent = variableArray{ii};

    fid = fopen(filePath);
    tline = fgetl(fid);
    tlineCounter = 1;

    % counter to keep track of repeated variable names
    varCtr = 1;

    % extract lines containing variables
    varLines = {};

    % loop line-by-line to find (multiple) instances of variable name
    while ischar(tline)
        if contains(tline, [varCurrent ' = '])
            if ~sum(tline(1:numel(varCurrent)) == varCurrent) == 0 % eliminate strings contained within other strings
                varLines{varCtr,1} = tline;
                varIdx(varCtr,1) = tlineCounter;
                varCtr = varCtr + 1;
            end
        end
        % Read next line in fid
        tline = fgetl(fid);
        tlineCounter = tlineCounter + 1;
    end
    fclose(fid);

    % find first instance of variable after object name
    varFirst = find(varIdx >= objectLine,1,'first');
    str = varLines{varFirst};

    % find values after: ' = '
    eqIdx = strfind(str, '=');
    str(1:eqIdx+1) = []; % +1 because of space after '='
    str(strfind(str, ',')) = [];
    
    % if variable is number(s):
    values{ii,1} = sscanf(str,'%f');

    % if variable value is string:
    if isempty(values{ii,1})
        values{ii,1} = sscanf(str,'%s');
    end
end


end  % get_goalc_parameters(...)