function fcmr_4dflow_postprocessing( reconDir, varargin )
%FCMR_4DFLOW_POSTPROCESSING  post-processing of 4D flow real-time data
%
%   FCMR_4DFLOW_POSTPROCESSING( reconDir, 'param', val )
%       Post-processing of 4D flow cine volumes and .vtk generation for 4D
%       visualisation in Paraview
%       - Loads 4D velocity component volumes and combines to single vector
%       - Applies background velocity drift correction
%       - Applies masks, i.e: blood pool masks
%       - Generates .vtk files for 4D visualisation in Paraview
%       - Generates .nii.gz files for in-plane visualisation in MRtrix
% 
%   Requires:
%       
%   Input:
%       reconDir            - str - local directory of fetal subjects
%
%   Optional Parameter-value Pairs:
%       velDir              directory containing SVR 4D velocity recon
%       fileExt             useful for running multiple postprocessing
%       useVelDriftCorr     velocity vector drift correction
%       bloodpoolMask       mask for cine volume
%       velMasks            cell array of masks for velocity volume -
%                           permits masking specific vessels only in
%                           4D flow volume and not the magnitude cine
%                           volume
%
%   Example usage:
%       fcmr_4dflow_postprocessing( 'C:\Users\tr17\Documents\Projects\PC_Fetal_CMR\Data\fcmr194' )
%
%   TODO:
%       - fix necessity for reslice_nii. Slight difference with SVRTK which means
%       reslicing of .nii necessary for MATLAB. Annoying having similarly
%       semi-duplicate files
%       - shift automatic masking from drift correction script
%
%   See also:
%       

% Tom Roberts (t.roberts@kcl.ac.uk)


%% Parse Input

default.velDir            = 'vel_vol';
default.fcmrNum           = [];
default.fileExt           = '';
default.useVelDriftCorr   = false;
default.bloodpoolMask     = 'mask_blood_pool';
default.velMasks          = {'mask_blood_pool'};

p = inputParser;

if  verLessThan('matlab','8.2')
    add_param_fn = @( parseobj, argname, defaultval, validator ) addParamValue( parseobj, argname, defaultval, validator );
else
    add_param_fn = @( parseobj, argname, defaultval, validator ) addParameter( parseobj, argname, defaultval, validator );
end

addRequired(  p, 'studyDir' );

add_param_fn( p, 'velDir', default.velDir, ...
        @(x) validateattributes( x, {'char'}, ...
        {}, mfilename ) );
add_param_fn( p, 'fcmrNum', default.fcmrNum, ...
        @(x) validateattributes( x, {'double'}, ...
        {}, mfilename ) );
add_param_fn( p, 'fileExt', default.fileExt, ...
        @(x) validateattributes( x, {'char'}, ...
        {}, mfilename ) );
add_param_fn( p, 'useVelDriftCorr', default.useVelDriftCorr, ...
        @(x) validateattributes( x, {'logical'}, ...
        {}, mfilename ) );
add_param_fn( p, 'bloodpoolMask', default.bloodpoolMask, ...
        @(x) validateattributes( x, {'char'}, ...
        {}, mfilename ) );
add_param_fn( p, 'velMasks', default.velMasks, ...
        @(x) validateattributes( x, {'cell'}, ...
        {}, mfilename ) );

parse( p, reconDir, varargin{:} );

velDir             = p.Results.velDir;
fcmrNum            = p.Results.fcmrNum;
fileExt            = p.Results.fileExt;
useVelDriftCorr    = p.Results.useVelDriftCorr;
bloodpoolMask      = p.Results.bloodpoolMask;
velMasks           = p.Results.velMasks;


%% Velocity volume polynomial correction
% NB: this automatically makes blood pool mask used further on
% TODO: make automatic masking separate from drift correction
if useVelDriftCorr == true
    fcmr_4dflow_drift_correction( reconDir , velDir );
end


%% Load cine volume
cd(reconDir);
cd('cine_vol');

if ~isfile('cine_vol-RESLICE.nii.gz')
    reslice_nii('cine_vol.nii.gz', 'cine_vol-RESLICE.nii.gz');
end
cine_nii = load_nii('cine_vol-RESLICE.nii.gz');


%% Number of frames in cine volume
nFrame = size(cine_nii.img,4);


%% Load velocity component volumes
cd(reconDir);
cd(velDir);

if strcmp(fileExt,'')
    velStr = 'velocity-final';
else
	velStr = ['velocity-final-' fileExt];
end

% reslice if not already performed
if ~isfile( [velStr '-RESLICE-0.nii.gz'] ) || ...
   ~isfile( [velStr '-RESLICE-1.nii.gz'] ) || ...
   ~isfile( [velStr '-RESLICE-2.nii.gz'] ) 
    reslice_nii([velStr '-0.nii.gz'], [velStr '-RESLICE-0.nii.gz']);
    reslice_nii([velStr '-1.nii.gz'], [velStr '-RESLICE-1.nii.gz']);
    reslice_nii([velStr '-2.nii.gz'], [velStr '-RESLICE-2.nii.gz']);
end

% load component vols
velx_nii = load_nii([velStr '-RESLICE-0.nii.gz']);
vely_nii = load_nii([velStr '-RESLICE-1.nii.gz']);
velz_nii = load_nii([velStr '-RESLICE-2.nii.gz']);

Vx = velx_nii.img;
Vy = vely_nii.img;
Vz = velz_nii.img;

% cm/s
Vx = 1e2 .* Vx;
Vy = 1e2 .* Vy; 
Vz = 1e2 .* Vz;


%% Load masks
cd(reconDir);
cd mask

%% Load magnitude cine blood pool mask
% NB: automatically created in fcmr_pc_velocity_correction.m

% reslice if necessary
if ~isfile([bloodpoolMask '-RESLICE.nii.gz'])
    reslice_nii([bloodpoolMask '.nii.gz'], [bloodpoolMask '-RESLICE.nii.gz']);
end

mask_bp = load_nii([bloodpoolMask '-RESLICE.nii.gz']);
mask_bp.img = double(mask_bp.img);
for ii = 1:nFrame; mask_bp.img(:,:,:,ii) = mask_bp.img(:,:,:,1); end


%% Load custom masks for 4D flow volume
% -default: blood pool mask, i.e: same as 4D cine volume
for mm = 1:numel(velMasks)      
    
    % current mask
    maskFileName = velMasks{mm};
    
    % open mask / reslice if necessary
    if ~isfile([maskFileName '-RESLICE.nii.gz'])
        reslice_nii([maskFileName '.nii.gz'], [maskFileName '-RESLICE.nii.gz']);
    end

    mask = load_nii([maskFileName '-RESLICE.nii.gz']);
    mask.img = double(mask.img);
    
    % fix mask_aorta/mask_IVC_SVC for fcmr194 (to do with extra voxels when re-slicing):
    fcmrNum = 194;
    if any(fcmrNum)
        if fcmrNum == 194 && strcmp(maskFileName,'mask_aorta') == 1 || fcmrNum == 194 && strcmp(maskFileName,'mask_IVC_SVC') == 1
            for tt = 1:nFrame; mask_re.img(:,:,:,tt) = imresize3(mask.img(:,:,:,tt),size(cine_nii.img(:,:,:,tt))); end
            mask.img = mask_re.img; clear mask_re;
        end
    end
    
    % initalise maskCombined
    if mm == 1
        maskCombined.img = zeros(size(mask.img));
    end

    maskCombined.img = maskCombined.img + mask.img;
%     maskCombined.img(:,:,1:64,:) = maskCombined.img(:,:,1:64,:) + mask.img(:,:,1:64,:);
    
end

mask.img = single(logical(maskCombined.img));
for ii = 1:nFrame; mask.img(:,:,:,ii) = mask.img(:,:,:,1); end


%% Resize
% FIXME: for some reason mask is 1 voxel larger than original volume...? 
for tt = 1:nFrame
    Vx_re(:,:,:,tt) = imresize3(Vx(:,:,:,tt),size(cine_nii.img(:,:,:,tt)));
    Vy_re(:,:,:,tt) = imresize3(Vy(:,:,:,tt),size(cine_nii.img(:,:,:,tt)));
    Vz_re(:,:,:,tt) = imresize3(Vz(:,:,:,tt),size(cine_nii.img(:,:,:,tt)));
    mask_re.img(:,:,:,tt) = imresize3(mask.img(:,:,:,tt),size(cine_nii.img(:,:,:,tt)));
    mask_bp_re.img(:,:,:,tt) = imresize3(mask_bp.img(:,:,:,tt),size(cine_nii.img(:,:,:,tt)));
end

Vx = Vx_re; clear Vx_re;
Vy = Vy_re; clear Vy_re;
Vz = Vz_re; clear Vz_re;
mask.img = mask_re.img; clear mask_re;
mask_bp.img = mask_bp_re.img; clear mask_bp_re;


%% Apply blood pool mask to cine_nii / custom masks to velocity volumes

% magnitude come
cine_masked_nii.img = double(cine_nii.img) .* mask_bp.img;

% velocity cine
Vx_masked = Vx .* mask.img;
Vy_masked = Vy .* mask.img;
Vz_masked = Vz .* mask.img;

Vmag_masked = sqrt(Vx_masked.^2 + Vy_masked.^2 + Vz_masked.^2);


%% Make meshgrids for VTK
[x,y,z] = meshgrid(1:size(Vx_masked,2),1:size(Vx_masked,1),1:size(Vx_masked,3));


%% Eliminate spurious values from Vmag
% improves Paraview visualisation
% hist(nonzeros(Vmag_masked(:)))
errIdx = find(Vmag_masked(:) > 100);
Vx_masked(errIdx) = 0; Vy_masked(errIdx) = 0; Vz_masked(errIdx) = 0; Vmag_masked(errIdx) = 0;


%% Set Low Velocity Threshold
% - cleans up volumes for viewing / fewer vectors to render
lowerBoundIdx = find(Vmag_masked(:) < 1);
Vx_masked(lowerBoundIdx) = 0; Vy_masked(lowerBoundIdx) = 0; Vz_masked(lowerBoundIdx) = 0; Vmag_masked(lowerBoundIdx) = 0;


%% Folder management if multiple velocity masks
%- string management if not just using mask_blood_pool
if sum(strcmp(velMasks,'mask_blood_pool')) == 1
    foldnameAppend = '';
else
    maskNames = velMasks;
    for ff = 1:numel(velMasks)
        maskNames{ff} = maskNames{ff}(6:end); 
    end
    foldnameAppend = ['_' strjoin(maskNames,'_')];
end

%% Write .vtk files for Paraview
cd(reconDir);
mkdir([velDir '_4d']);
cd([velDir '_4d']);

if strcmp(fileExt,'')
    mkdir(['paraview' foldnameAppend]);
    cd(['paraview' foldnameAppend])
else
	mkdir(['paraview_' fileExt foldnameAppend]);
    cd(['paraview_' fileExt foldnameAppend]);
end


disp('Writing .vtk files ... ');

% make .vtk files
for tt = 1:nFrame
    
    % magnitude cine volume
    vtkwrite(['cine_vol_masked_t' num2str(tt-1) '.vtk'], ...
             'structured_grid', y, x, z, ... 
             'scalars', 'magnitude_intensity', cine_masked_nii.img(:,:,:,tt) );
         
    % velocity cine volume
    vtkwrite(['vel_vol_masked_VxVy-Vz_t' num2str(tt-1) '.vtk'], ...
             'structured_grid', y, x, z, ... 
             'vectors', 'vector_field', Vx_masked(:,:,:,tt), Vy_masked(:,:,:,tt), -1 .* Vz_masked(:,:,:,tt), ...
             'scalars', 'velocity_magnitude', Vmag_masked(:,:,:,tt) );  
         
%     %%% Correct configuration of VxVyVz in Flow Phantom 6 
%     vtkwrite(['vel_vol_masked_VyVx-Vz_t' num2str(tt-1) '.vtk'], ...
%              'structured_grid', x, y, z, ... 
%              'vectors', 'vector_field', Vy_masked(:,:,:,tt), Vx_masked(:,:,:,tt), -1 .* Vz_masked(:,:,:,tt), ...
%              'scalars', 'velocity_magnitude', Vmag_masked(:,:,:,tt) );  

end

disp('Finished making .vtk files ... ');


%% Write .nii for MRtrix
cd(reconDir);
cd([velDir '_4d']);

if strcmp(fileExt,'')
    mkdir(['mrtrix' foldnameAppend]);
    cd(['mrtrix' foldnameAppend])
    
    % get .nii to use as basis
    Vx3D_nii = load_untouch_nii([reconDir '\' velDir '\velocity-final-RESLICE-0.nii.gz']);
else
	mkdir(['mrtrix_' fileExt foldnameAppend]);
    cd(['mrtrix_' fileExt foldnameAppend]);
    
    % get .nii to use as basis
    Vx3D_nii = load_untouch_nii([reconDir '\' velDir '\velocity-final-' fileExt '-RESLICE-0.nii.gz']);
end


disp('Writing .nii.gz files ... ');

for tt = 1:nFrame
    v3D = Vx3D_nii;
    v3D.img = cat(4,Vx_masked(:,:,:,tt), Vy_masked(:,:,:,tt), -1 .* Vz_masked(:,:,:,tt));
    
    v3D.hdr.dime.dim(1) = 4;
    v3D.hdr.dime.dim(2:5) = [size(Vx_masked,1) size(Vx_masked,2) size(Vx_masked,3) 3];
    v3D.hdr.dime.pixdim(5) = 0;
    
    save_untouch_nii(v3D,['VxVy-Vz_t' num2str(tt-1) '.nii.gz']);
end

disp('Finished making .nii.gz files ...');


%% Return to reconDir
cd(reconDir);

disp('Post-processing complete ...');


end %fcmr_4dflow_postprocessing(...)
