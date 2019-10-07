function values = get_pcmr_orientation_parameters( gcFilePath )

%GET_PCMR_ORIENTATION_PARAMETERS extract important PC-bSSFP parameters from goalc.txt file
%
%   GET_PCMR_ORIENTATION_PARAMETERS( goalcFile, objectStr, variableArray )
% 
%   Output:
%       values              - structure containing specified variables and
%                             associated values:
%                             O_MATRIX/GR(readout)/GR(slice-select)
%
%   Input:
%       gcFilePath          - file path to input goalc.txt file
%
%   Example:
%       values = get_pcmr_orientation_parameters( <../goalc.txt> )
%
%   Usage:
%       Use with Josh's mrecon function: mrecon_writegoalc2txt.m
%       Use with: get_goalc_parameters.m

% TAR (t.roberts@kcl.ac.uk)


%% Initialise
values = struct;


%% O_MATRIX - contains slice orientation MPS->xyz

object = 'O_MATRIX `locations[0]';
varArray = { 'm_orient', 'p_orient', 's_orient' }; % variables to collect from object

for ii = 1:numel(varArray)       
    values.(varArray{ii}) = cell2mat(get_goalc_parameters( gcFilePath, object, varArray{ii} ))';
end


%% Variables to collect from each GR object
varArray = { 'str', 'flow1', 'flow2', 'ori_x', 'ori_y', 'ori_z' }; 


%% Measurement gradient objects

% Collect variables from objects
objectArray = { 'GR `mc[0]', 'GR `m[0]', 'GR `m[3]' };
structNames = { 'mc0', 'm0', 'm3' }; %strings in objectArray incompatiable with MATLAB structures

for oo = 1:numel(objectArray)
    
    % prepare structure names for current object
    structNameArray = { [structNames{oo} '_str'],  ...
                        [structNames{oo} '_flow1'], ...
                        [structNames{oo} '_flow2'], ...
                        [structNames{oo} '_ori_x'], ...
                        [structNames{oo} '_ori_y'], ...
                        [structNames{oo} '_ori_z'] }; % variables to collect from object

	% get variables for current object
    for ii = 1:numel(varArray)            
            values.(structNameArray{ii}) = cell2mat(get_goalc_parameters( gcFilePath, objectArray{oo}, varArray{ii} ));
    end
end


%% Slice-select gradient objects

% Collect variables from objects
objectArray = { 'GR `s_ex', 'GR `r[0]' };
structNames = { 's_ex', 'r0' };

for oo = 1:numel(objectArray)
    
    % prepare structure names for current object
    structNameArray = { [structNames{oo} '_str'],  ...
                        [structNames{oo} '_flow1'], ...
                        [structNames{oo} '_flow2'], ...
                        [structNames{oo} '_ori_x'], ...
                        [structNames{oo} '_ori_y'], ...
                        [structNames{oo} '_ori_z'] }; % variables to collect from object

	% get variables for current object
    for ii = 1:numel(varArray)            
            values.(structNameArray{ii}) = cell2mat(get_goalc_parameters( gcFilePath, objectArray{oo}, varArray{ii} ));
    end
end


end % get_pcmr_orientation_parameters(...)