function D = summarise_transformations( recondir, resultstable, varargin )
%SUMMARISE_TRANSFORMATIONS  summarises volumetric reconstruction
%transformation results
%
%   D = SUMMARISE_TRANSFORMATIONS( recondir, resultstable ) returns results
%   in structure D given directory for volumetric reconstruction data and
%   table of results in .tsv file produced by read_info_tsv().
% 
%   SUMMARISE_TRANSFORMATIONS( ..., 'name', value ) specifies optional 
%   input arguments. See code for name-value pairs.
% 
%   See also read_info_tsv.

%   jfpva (joshua.vanamerom@kcl.ac.uk) 


%% Optional Input Argument Default Values

default.niftiPath       = fullfile(fileparts(mfilename('fullpath')),'..','lib','nifti');  % path to NIfTI toolbox
default.irtkPath        = '/usr/local/build/bin/';  % path to IRTK bin 


%% Parse Input

p = inputParser;

if  verLessThan('matlab','8.2')
    add_param_fn = @( parseobj, argname, defaultval, validator ) addParamValue( parseobj, argname, defaultval, validator );
else
    add_param_fn = @( parseobj, argname, defaultval, validator ) addParameter( parseobj, argname, defaultval, validator );
end

addRequired(   p, 'recondir', ... 
    @(x) validateattributes( x, {'char'}, {'vector'}, mfilename) );

addRequired(   p, 'resultstable', ... 
    @(x) validateattributes( x, {'table'}, {}, mfilename) );

add_param_fn(   p, 'niftipath', default.niftiPath, ...
    @(x) validateattributes( x, {'char'}, {'vector'}, mfilename) );

add_param_fn(   p, 'irtkpath', default.irtkPath, ...
    @(x) validateattributes( x, {'char'}, {'vector'}, mfilename) );

parse( p, recondir, resultstable, varargin{:} );

reconDir            = p.Results.recondir;
T                   = p.Results.resultstable;
niftiPath           = p.Results.niftipath;
irtkPath            = p.Results.irtkpath;


%% Dependencies

% Temporarily add Matlab toolboxes
origPath  = path;
resetPath = onCleanup( @() path( origPath ) );
addpath( genpath( niftiPath ) )

% Temporarily add IRTK to matlab environment
origEnv  = getenv('PATH');
resetEnv = onCleanup( @() setenv( 'PATH', origEnv ) );
setenv( 'PATH', [ origEnv ':' irtkPath ] );  


%% Read Stack to World Transformations and Stack Image Coordinates

S = dir( fullfile( reconDir, 'stack*.nii.gz' ) );

stackCoord      = cell( 1, size(T,1) );

Astack2world    = cell( 1, numel(S) );

sliceIndex      = nan( size( stackCoord ) );
stackIndex      = nan( size( stackCoord ) );

for iSt = 1:numel(S)
    
    % identify stack file
    stackFilePath = fullfile( S(iSt).folder, S(iSt).name );

    % load stack images
	N = load_untouch_nii( stackFilePath );

    % find image voxel coordinates
    ind = find(N.img~=-1); 
    [i,j,k,d] = ind2sub( size(N.img), ind );

    % adjust image coordinates to be zero-indexed
    i = i-1; 
    j = j-1; 
    k = k-1;
    d = d-1;
    
    % assemble image coordinates as 4x1 vector for each voxel, 
    % with 1 in fourth position to allow multiplication by augmented 
    % transformation matrices
    inputIndex    = T.InputIndex(T.StackIndex==(iSt-1)); 
    stackLocIndex = T.StackLocIndex(T.StackIndex==(iSt-1)); 
    stackDynIndex = T.StackDynIndex(T.StackIndex==(iSt-1)); 
    locIndex      = T.LocIndex(T.StackIndex==(iSt-1)); 
    for iStK = 1:numel(inputIndex)
        ind = ( k==stackLocIndex(iStK) & d==stackDynIndex(iStK) );
        stackCoord{inputIndex(iStK)+1} = [ i(ind).'; j(ind).'; k(ind).'; ones(1,sum(ind)) ];
        stackIndex(inputIndex(iStK)+1) = iSt;
        sliceIndex(inputIndex(iStK)+1) = locIndex(iStK)+1;
    end
    
    % extract stack to world transformation matrices
    cmd  = sprintf( '%s %s | grep -A 5 "Image to world matrix" | tail -4', fullfile( irtkPath, 'info' ), stackFilePath );
    [ status, result ] = system( cmd );
    if ( (status==0) && ~isempty(result) )
        try
            Astack2world{iSt} = eval( strcat( '[', result, ']' ) );
        catch
            fprintf( 'Could not read %s image-to-world transformation\n', stackFilePath )
            Astack2world{iSt} = nan(4);
        end
    else
        fprintf( 'Could not read %s image-to-world transformation\n', stackFilePath )
        Astack2world{iSt} = nan(4);
    end
      
end


%% Read Estimated Image-Volume Transformations

Aim2vol = cell( 1, size(T,1) );

for iK = 1:size(T,1)
    dofFilePath = fullfile( reconDir, 'transformations', sprintf( 'transformation%05i.dof', iK-1 ) );
	cmd  = sprintf( '%s %s | tail -4', fullfile( irtkPath, 'dof2mat' ), dofFilePath );
    [ status, result ] = system( cmd );
    if ( (status==0) && ~isempty(result) )
        try
            Aim2vol{iK} = eval( strcat( '[', result, ']' ) );
        catch
            fprintf( 'Could not read %s\n', dofFilePath )
            Aim2vol{iK} = nan(4);
        end
    else
        fprintf( 'Could not read %s\n', dofFilePath )
        Aim2vol{iK} = nan(4);
    end  
end


%% Calculate World Coordinates

worldCoord = cell( 1, size(T,1) );

for iK = 1:size(T,1)
    worldCoord{iK} = Astack2world{stackIndex(iK)} * stackCoord{iK};
end


%% Calculate Slice Transformations

Aslice2vol = cell( 1, numel( unique( T.LocIndex ) ) );

for iSl = unique(sliceIndex)
    iK = find(sliceIndex==iSl);
    A = Aim2vol(iK);
    wgt = T.Weight(iK);  
    while( sum(wgt) == 0 )
        wgtOffset = 1e-13;
        wgt = wgt + wgtOffset;  % add small value to weights in case that all pFrame for slice are zero 
        warning( 'Sum of weights (pframe for slice l) <= 0 for l = %i; adding %g to weights for all frames in slice', iSl-1, wgtOffset )
    end
    Aslice2vol{iSl} = calc_frechet_mean( A, wgt );
end


%% Calculate Displacements and Deviations

calc_dist = @(A1,A2,worldCoordinates) sqrt( sum( ( A1*worldCoordinates-A2*worldCoordinates ).^2, 1 ) );
calc_disp = @(dist) mean( dist );
calc_dev  = @(dist) mean( dist );

% Image Frame Displacement
dispIm = nan( 1, size(T,1) );
distDispIm = cell( 1, size(T,1) );
for iK = 1:size(T,1)
    distDispIm{iK} = calc_dist( eye(4), Aim2vol{iK}, worldCoord{iK} );
    dispIm(iK)  = calc_disp( distDispIm{iK} );
end

% Global Displacement
dist = [];
for iK = 1:size(T,1)
    dist = cat( 2, dist, calc_dist( eye(4), Aim2vol{iK}, worldCoord{iK} ) ); 
end
dispAll = mean( dist );

% Image Frame Deviations
devIm = nan( 1, size(T,1) );
distDevIm = cell( 1, size(T,1) );
for iK = 1:size(T,1)
    distDevIm{iK} = calc_dist( Aslice2vol{sliceIndex(iK)}, Aim2vol{iK}, worldCoord{iK} );
    devIm(iK) = calc_dev( distDevIm{iK} );
end

% Slice Deviations
devSl = nan( 1, numel(unique(sliceIndex)) );
for iSl = unique(sliceIndex)
    d = [];
    for iK = find(sliceIndex==iSl)
        d = cat( 2, d, distDevIm{iK} );
    end
    devSl(iSl) = calc_dev( d );
end

% Global Deviation
d = [];
for iSl = unique(sliceIndex)
    for iK = find(sliceIndex==iSl)
        d = cat( 2, d, distDevIm{iK} );
    end
end
devAll = calc_dev( d );


%% Combine Results as Output Structure

D.stackCoord        = stackCoord;
D.worldCoord        = worldCoord;
D.sliceIndex        = sliceIndex;
D.stackIndex        = stackIndex;
D.Astack2world      = Astack2world;
D.Aim2vol           = Aim2vol;
D.Aslice2vol        = Aslice2vol;
D.distDispIm        = distDispIm;
D.dispIm            = dispIm;
D.dispAll           = dispAll;
D.distDevIm         = distDevIm;
D.devIm             = devIm;
D.devSl             = devSl;
D.devAll            = devAll;


end  % summarise_transformations(...)