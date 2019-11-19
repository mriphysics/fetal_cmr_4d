function get_cine_vol_force_exclude(reconDir)

%GET_CINE_VOL_FORCE_EXCLUDE  get excluded frame numbers following cine_vol recon
%
% Usage: must have cine_vol folder reconstructed. Reads excluded frames from
% log-evaluation.txt and writes to new text file for use with
% recon_vel_vol.bash
%
% Example:
%   get_cine_vol_force_exclude('~/path/to/reconDir/')
%
%    See also fcmr_4dflow_preprocessing.m.

%   tar (t.roberts@kcl.ac.uk) 

cd(reconDir);    
cd cine_vol

% Open cine_vol log-evaluation.txt
fid = fopen('log-evaluation.txt');
C = textscan(fid,'%s');
C = C{1,1};
fclose(fid);

excludedIdx = find(strcmp(C,'Excluded'),1,'last');
outsideIdx  = find(strcmp(C,'Outside'),1,'last');
totalIdx    = find(strcmp(C,'Total:'),2,'last');

% Get Excluded Slices
if ~isempty(excludedIdx+2:totalIdx(1)-1)

    ctr = 1;
    for ii = excludedIdx+2:totalIdx(1)-1
        excludedVals(ctr) = str2double(cell2mat(C(ii)));
        ctr = ctr + 1;
    end

    if ( numel(excludedVals) ~= str2double(cell2mat(C(totalIdx(1)+1))) )
        error('Number of Excluded Slices does not match.');
    end
    
elseif isempty(excludedIdx+2:totalIdx(1)-1)
    
    excludedVals = [];
    
end
    
% Get Outside Slices
if ~isempty(outsideIdx+2:totalIdx(2)-1)

    ctr = 1;
    for ii = outsideIdx+2:totalIdx(2)-1
        outsideVals(ctr) = str2double(cell2mat(C(ii)));
        ctr = ctr + 1;
    end

    if ( numel(outsideVals) ~= str2double(cell2mat(C(totalIdx(2)+1))) )
        error('Number of Outside Slices does not match.');
    end
    
elseif isempty(outsideIdx+2:totalIdx(2)-1)
    
    outsideVals = [];
    
end

% Total Combined Number of Outisde/Excluded Slices
totalVals = numel(excludedVals) + numel(outsideVals);

% Write Excluded cine_vol Slices To Text File
cd ../data
fid = fopen( 'force_exclude_cine_vol.txt' , 'w' );
% fprintf( fid, '%d ', totalVals ); % Not required as calculated in Shell script
fprintf( fid, '%d ', sort([excludedVals,outsideVals]) );
fclose( fid );

disp(['force_exclude_cine_vol.txt written to: ' fullfile( reconDir, 'data' ) ]);

% fn end
end