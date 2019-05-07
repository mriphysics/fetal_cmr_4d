function S = preproc( reconDir )
%PREPROC  preprocess data for 4D reconstruction.
%
%   S = PREPROC( reconDir ) reads files from reconDir and returns data structure S.

%   JFPvA (joshua.vanamerom@kcl.ac.uk)


%% Init

dataDir     = fullfile( reconDir, 'data' );
maskDir     = fullfile( reconDir, 'mask' );
ktreconDir  = fullfile( reconDir, 'ktrecon' );
hrRange     = [105,180];
nFrameCine  = 25;
delta       = 150;


%% Identify Data and Get Parameters

% Identify Dynamic MR Image Series
rltFileList       = dir( fullfile( dataDir, '*_rlt_ab.nii.gz' ) );

% Get Number of Stacks
nStack = numel(rltFileList);

% Initialise Stack Struct
S = struct([]);

% Read Data Into Stack Struct
for iStk = 1:nStack

    % Identify Files
    S(iStk).desc          = strrep( rltFileList(iStk).name, '_rlt_ab.nii.gz', '' );
    S(iStk).rltAbFile     = fullfile( rltFileList(iStk).folder,  rltFileList(iStk).name  );
    S(iStk).rltReFile     = fullfile( ktreconDir, strrep( rltFileList(iStk).name, 'ab', 're' ) );
    S(iStk).rltImFile     = fullfile( ktreconDir, strrep( rltFileList(iStk).name, 'ab', 'im' ) );
    S(iStk).rltMatFile    = fullfile( ktreconDir, sprintf( '%s_rlt_recon.mat', S(iStk).desc ) );
    S(iStk).rltParamFile  = fullfile( ktreconDir, sprintf( '%s_rlt_parameters.mat', S(iStk).desc ) );
    S(iStk).dcAbFile      = fullfile( dataDir, sprintf( '%s_dc_ab.nii.gz', S(iStk).desc ) );
    S(iStk).slwAbFile     = fullfile( ktreconDir, sprintf( '%s_slw_ab.nii.gz', S(iStk).desc ) );
    S(iStk).trnAbFile     = fullfile( ktreconDir, sprintf( '%s_trn_ab.nii.gz', S(iStk).desc ) );
    S(iStk).maskHeartFile = fullfile( maskDir, sprintf( '%s_mask_heart.nii.gz', S(iStk).desc ) );
       
    % Load Parameters
    M = matfile( S(iStk).rltParamFile );
    P = M.PARAM;

    % Extract Parameters
    S(iStk).nLoc             = P.Timing.numLoc;
    S(iStk).sliceThickness   = P.Scan.RecVoxelSize(3);
    
    % Load NIfTI
    R = load_untouch_nii( S(iStk).rltAbFile );
    S(iStk).niiHdr = R.hdr;

    % Separate slices
    for iLoc = 1:S(iStk).nLoc

        % Dynamic Image Series
        S(iStk).frameDuration   = P.Timing.frameDuration;
        S(iStk).tFrame{iLoc}    = P.Timing.sliceTime(iLoc) + S(iStk).frameDuration * (0:(P.Encoding.NrDyn(1)-1));
        
    end
    
end


%% Set Up Additional Data

filePath = fullfile( dataDir, 'tgt_stack_no.txt' );
if ~exist( filePath, 'file' )
  if ~exist( 'iStkTgt', 'var')
    iStkTgt = 1;
  end
    fid = fopen( filePath , 'w' );
    fprintf( fid, num2str(iStkTgt) );
    fclose( fid );
end

filePath = fullfile( dataDir, 'slice_thickness.txt' );
fid = fopen( filePath , 'w' );
fprintf( fid, '%g ', [S.sliceThickness] );
fclose( fid );

filePath = fullfile( dataDir, 'force_exclude_stack.txt' );
if ~exist( filePath, 'file' )
    fid = fopen( filePath , 'w' );
    fprintf( fid, '' );
    fclose( fid );
end

filePath = fullfile( dataDir, 'force_exclude_slice.txt' );
if ~exist( filePath, 'file' )
    fid = fopen( filePath , 'w' );
    fprintf( fid, '' );
    fclose( fid );
end

filePath = fullfile( dataDir, 'force_exclude_frame.txt' );
if ~exist( filePath, 'file' )
    fid = fopen( filePath , 'w' );
    fprintf( fid, '' );
    fclose( fid );
end


end  % preproc(...)

