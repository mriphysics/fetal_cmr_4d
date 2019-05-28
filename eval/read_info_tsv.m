function T = read_info_tsv( infoFilename, logFilename, startRow, endRow)
%READ_INFO_TSV  import numeric data from a text file as a matrix
%
%   T = READ_INFO_TSV(INFOFILENAME) Reads data from text file 
%   INFOFILENAME for the default selection.
%
%   T = READ_INFO_TSV(INFOFILENAME,LOGFILENAME) Reads additional data from 
%   log file LOGFILENAME.
%
%   T = READ_INFO_TSV(FILENAME, LOGFILENAME, STARTROW, ENDROW) Reads data from
%   rows STARTROW through ENDROW of text file FILENAME.
%
% Example:
%   T = read_info_tsv('info.tsv','log-reconstruction.txt');
%
%    See also TEXTSCAN.

%   jfpva (joshua.vanamerom@kcl.ac.uk) 


%% Initialize variables.
delimiter = '\t';
if nargin<=2
    startRow = 2;
    endRow = inf;
end

%% Read columns of data as text:
% For more information, see the TEXTSCAN documentation.
formatSpec = '%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%[^\n\r]';

%% Open the text file.
fileID = fopen(infoFilename,'r');

%% Read columns of data according to the format.
% This call is based on the structure of the file used to generate this
% code. If an error occurs for a different file, try regenerating the code
% from the Import Tool.
dataArray = textscan(fileID, formatSpec, endRow(1)-startRow(1)+1, 'Delimiter', delimiter, 'TextType', 'string', 'HeaderLines', startRow(1)-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
for block=2:length(startRow)
    frewind(fileID);
    dataArrayBlock = textscan(fileID, formatSpec, endRow(block)-startRow(block)+1, 'Delimiter', delimiter, 'TextType', 'string', 'HeaderLines', startRow(block)-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
    for col=1:length(dataArray)
        dataArray{col} = [dataArray{col};dataArrayBlock{col}];
    end
end

%% Close the text file.
fclose(fileID);

%% Convert the contents of columns containing numeric text to numbers.
% Replace non-numeric text with NaN.
raw = repmat({''},length(dataArray{1}),length(dataArray)-1);
for col=1:length(dataArray)-1
    raw(1:length(dataArray{col}),col) = mat2cell(dataArray{col}, ones(length(dataArray{col}), 1));
end
numericData = NaN(size(dataArray{1},1),size(dataArray,2));

for col=[1,2,3,4,5,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28]
    % Converts text in the input cell array to numbers. Replaced non-numeric
    % text with NaN.
    rawData = dataArray{col};
    for row=1:size(rawData, 1)
        % Create a regular expression to detect and remove non-numeric prefixes and
        % suffixes.
        regexstr = '(?<prefix>.*?)(?<numbers>([-]*(\d+[\,]*)+[\.]{0,1}\d*[eEdD]{0,1}[-+]*\d*[i]{0,1})|([-]*(\d+[\,]*)*[\.]{1,1}\d+[eEdD]{0,1}[-+]*\d*[i]{0,1}))(?<suffix>.*)';
        try
            result = regexp(rawData(row), regexstr, 'names');
            numbers = result.numbers;
            
            % Detected commas in non-thousand locations.
            invalidThousandsSeparator = false;
            if numbers.contains(',')
                thousandsRegExp = '^\d+?(\,\d{3})*\.{0,1}\d*$';
                if isempty(regexp(numbers, thousandsRegExp, 'once'))
                    numbers = NaN;
                    invalidThousandsSeparator = true;
                end
            end
            % Convert numeric text to numbers.
            if ~invalidThousandsSeparator
                numbers = textscan(char(strrep(numbers, ',', '')), '%f');
                numericData(row, col) = numbers{1};
                raw{row, col} = numbers{1};
            end
        catch
            raw{row, col} = rawData{row};
        end
    end
end


%% Split data into numeric and string columns.
rawNumericColumns = raw(:, [1,2,3,4,5,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28]);
rawStringColumns = string(raw(:, 6));


%% Replace non-numeric cells with NaN
R = cellfun(@(x) ~isnumeric(x) && ~islogical(x),rawNumericColumns); % Find non-numeric cells
rawNumericColumns(R) = {NaN}; % Replace non-numeric cells

%% Create output variable
T = table;
T.StackIndex = cell2mat(rawNumericColumns(:, 1));
T.StackLocIndex = cell2mat(rawNumericColumns(:, 2));
T.StackDynIndex = cell2mat(rawNumericColumns(:, 3));
T.LocIndex = cell2mat(rawNumericColumns(:, 4));
T.InputIndex = cell2mat(rawNumericColumns(:, 5));
T.File = rawStringColumns(:, 1);
T.Scale = cell2mat(rawNumericColumns(:, 6));
T.StackFactor = cell2mat(rawNumericColumns(:, 7));
T.Time = cell2mat(rawNumericColumns(:, 8));
T.TemporalResolution = cell2mat(rawNumericColumns(:, 9));
T.CardiacPhase = cell2mat(rawNumericColumns(:, 10));
T.ReconCardPhaseIndex = cell2mat(rawNumericColumns(:, 11));
T.Included = cell2mat(rawNumericColumns(:, 12));
T.Excluded = cell2mat(rawNumericColumns(:, 13));
T.Outside = cell2mat(rawNumericColumns(:, 14));
T.Weight = cell2mat(rawNumericColumns(:, 15));
T.MeanDisplacement = cell2mat(rawNumericColumns(:, 16));
T.MeanDisplacementX = cell2mat(rawNumericColumns(:, 17));
T.MeanDisplacementY = cell2mat(rawNumericColumns(:, 18));
T.MeanDisplacementZ = cell2mat(rawNumericColumns(:, 19));
T.WeightedMeanDisplacement = cell2mat(rawNumericColumns(:, 20));
T.TRE = cell2mat(rawNumericColumns(:, 21));
T.TranslationX = cell2mat(rawNumericColumns(:, 22));
T.TranslationY = cell2mat(rawNumericColumns(:, 23));
T.TranslationZ = cell2mat(rawNumericColumns(:, 24));
T.RotationX = cell2mat(rawNumericColumns(:, 25));
T.RotationY = cell2mat(rawNumericColumns(:, 26));
T.RotationZ = cell2mat(rawNumericColumns(:, 27));

%% Replace invalid values

T.MeanDisplacement(T.MeanDisplacement==-1) = NaN;
T.WeightedMeanDisplacement(T.WeightedMeanDisplacement==-1) = NaN;
T.TRE(T.TRE==-1) = NaN;

%% Read Small Image Frames from Log File

if exist( 'logFilename', 'var' )
    cmd = sprintf( 'grep "Small slices:" %s | tail -1', logFilename );
    [ ~, result ] = system( cmd );
    smallImages = str2num( erase( erase( result, newline ), 'Small slices: ' ) ); %#ok<*ST2NM>
    T.Small = zeros( size( T.Outside ) );
    T.Small(smallImages+1) = 1;
end

%% Read Image Potentials from Log File

if exist( 'logFilename', 'var' )
    cmd = sprintf( 'grep "Slice potentials: " %s | tail -1', logFilename );
    [ ~, result ] = system( cmd );
    if ~isempty(result)
        potential = str2num( erase( erase( result, newline ), 'Slice potentials: ' ) );
        T.Potential = potential(:);
    else
        T.Potential = nan( size(T,1), 1 );
    end
end

end  % read_info_tsv(...)

