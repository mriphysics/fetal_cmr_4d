function fcmr_4dflow_drift_correction( reconDir , velDir, bpThreshold )

%% fcmr_4dflow_drift_correction()
% Apply polynomial to velocity volumes
%
% Tom Roberts (t.roberts@kcl.ac.uk)

%% Admin

% Set default bloodpool threshold: higher = less blood pool
if nargin < 3
    bpThreshold = 100;
elseif nargin < 2
    bpThreshold = 100;
    velDir = 'vel_vol';
end

cd(reconDir);


%% Load data

% search for resliced .nii
cd cine_vol
% resliceFiles = dir('cine*RESLICE*nii.gz');

% if contains( resliceFiles.name , 'RESLICE' )
%     
%     warning('Using resliced cine_vol and velocity-final nifti files ... ');
%     
%     cineName = 'cine_vol-RESLICE';
%     cv = load_untouch_nii([cineName '.nii.gz']);
%     
%     cd ..
%     cd(velDir)
%     
%     velName = 'velocity-final-RESLICE';    
%     v0 = load_untouch_nii([velDir '/' velName '-0.nii.gz']);
%     v1 = load_untouch_nii([velDir '/' velName '-1.nii.gz']);
%     v2 = load_untouch_nii([velDir '/' velName '-2.nii.gz']);
% 
% else
    
    cineName = 'cine_vol';
    cv = load_untouch_nii([cineName '.nii.gz']);
    
    cd ..
    cd(velDir)
    
    velName = 'velocity-final';    
    v0 = load_untouch_nii([velDir '/' velName '-0.nii.gz']);
    v1 = load_untouch_nii([velDir '/' velName '-1.nii.gz']);
    v2 = load_untouch_nii([velDir '/' velName '-2.nii.gz']);
    
% end


%% Resample cine if dimensions do not exactly match velocity volume
if size(cv.img,1) ~= size(v0.img,1) || size(cv.img,2) ~= size(v0.img,2) || size(cv.img,3) ~= size(v0.img,3)
    
    warning('Resizing cine_vol to match velocity volumes ... ');
    
    for tt = 1:size(cv.img,4)
        cv_re.img(:,:,:,tt) = imresize3( cv.img(:,:,:,tt) , size(v0.img(:,:,:,1)) );
    end
    cv.img = cv_re.img;
end


%% Mask blood pool
% bpThreshold = 110; % 'blood pool' threshold - arbitrary number empirically determined
% BP = cv.img>bpThreshold;

BP = ~(cv.img<prctile(cv.img(:),90)); % percentile method

% View full cine_vol / blood pool only / inverse
montage_RR(cv.img(:,:,:,1),'gray',[0,150]); close;
montage_RR(BP(:,:,:,1) .* cv.img(:,:,:,1),'gray',[0,150]); close;
montage_RR(~BP(:,:,:,1) .* cv.img(:,:,:,1),'gray',[0,150]);


%% Mask background (and velocity 'rim') and in v0/v1/v2 --- currently set to -15
% BK = cv.img<0;
BK = v0.img==-15;
v0.img(BK) = 0; v1.img(BK) = 0; v2.img(BK) = 0;


%% Apply temporal mask
% - ONLY allow voxels which never see blood pool through time
NBP = logical(~BP & ~BK);           % non-blood pool
MT = sum(NBP,4);                    % sum masks through time
MT( MT < size(cv.img,4) ) = 0;      % set any voxels which see blood pool = 0
MT = logical(MT);

% Find average velocity through time in non-blood pool regions
v0_mean_time.img = mean( v0.img , 4  ) .* MT;
v1_mean_time.img = mean( v1.img , 4  ) .* MT;
v2_mean_time.img = mean( v2.img , 4  ) .* MT;

% Find median velocity through time in non-blood pool regions
v0_median_time.img = median( v0.img , 4 ) .* MT;
v1_median_time.img = median( v1.img , 4 ) .* MT;
v2_median_time.img = median( v2.img , 4 ) .* MT;

% STDEV velocity through time in non-blood pool regions
v0_std_time.img = std( v0.img , [], 4  ) .* MT;
v1_std_time.img = std( v1.img , [], 4  ) .* MT;
v2_std_time.img = std( v2.img , [], 4  ) .* MT;


%% Apply polynomial to DC velocity volume

disp('Applying polynomial correction to velocity volumes ...');

polyOrder = 3;
M = logical(MT);

% fit polynomial
[xdc,ydc,zdc] = meshgrid( 1:size(cv.img,2) , 1:size(cv.img,1) , 1:size(cv.img,3) );
model0 = polyfitn( [xdc(M(:)),ydc(M(:)),zdc(M(:))], v0_median_time.img(M(:)), polyOrder );
model1 = polyfitn( [xdc(M(:)),ydc(M(:)),zdc(M(:))], v1_median_time.img(M(:)), polyOrder );
model2 = polyfitn( [xdc(M(:)),ydc(M(:)),zdc(M(:))], v2_median_time.img(M(:)), polyOrder );

V0 = reshape( polyvaln( model0, [xdc(:),ydc(:),zdc(:)] ), [size(cv.img,1),size(cv.img,2),size(cv.img,3)] );
V1 = reshape( polyvaln( model1, [xdc(:),ydc(:),zdc(:)] ), [size(cv.img,1),size(cv.img,2),size(cv.img,3)] );
V2 = reshape( polyvaln( model2, [xdc(:),ydc(:),zdc(:)] ), [size(cv.img,1),size(cv.img,2),size(cv.img,3)] );

% Subtract polynomial from each timepoint
for tt = 1:size(v0.img,4)
    v0_corr.img(:,:,:,tt) = v0.img(:,:,:,tt) - V0;
    v1_corr.img(:,:,:,tt) = v1.img(:,:,:,tt) - V1;
    v2_corr.img(:,:,:,tt) = v2.img(:,:,:,tt) - V2;
end


%% Apply background mask
v0_corr.img = v0_corr.img .* ~BK;
v1_corr.img = v1_corr.img .* ~BK;
v2_corr.img = v2_corr.img .* ~BK;


%% Compare original velocity volumes to corrected volumes
imtar([v0.img(:,:,16,20) V0(:,:,16).*~BK(:,:,16) v0_corr.img(:,:,16,20); ...
       v0_median_time.img(:,:,16) v0_std_time.img(:,:,16) zeros(size(v0_std_time.img(:,:,16)))],-0.1,0.1);
close;

% View raw and polynomial corrected velocity volumes
implay_RR([v0.img(:,:,16,:) v0_corr.img(:,:,16,:) v0.img(:,:,16,:)-v0_corr.img(:,:,16,:); ...
           v1.img(:,:,16,:) v1_corr.img(:,:,16,:) v1.img(:,:,16,:)-v1_corr.img(:,:,16,:); ...
           v2.img(:,:,16,:) v2_corr.img(:,:,16,:) v2.img(:,:,16,:)-v2_corr.img(:,:,16,:); ...
           v0.img(:,:,16,:)+v1.img(:,:,16,:)+v2.img(:,:,16,:) v0_corr.img(:,:,16,:)+v1_corr.img(:,:,16,:)+v2_corr.img(:,:,16,:) ( v0.img(:,:,16,:)+v1.img(:,:,16,:)+v2.img(:,:,16,:) ) - ( v0_corr.img(:,:,16,:)+v1_corr.img(:,:,16,:)+v2_corr.img(:,:,16,:) )], ...
           'jet',[-0.2,0.2]);
       
       
%% Save polynomial corrected .nii
v0.img = v0_corr.img; v1.img = v1_corr.img; v2.img = v2_corr.img;

cd(reconDir);

% save_untouch_nii(v0,[velDir '/' velName '-bkCorr-0.nii.gz']);
% save_untouch_nii(v1,[velDir '/' velName '-bkCorr-1.nii.gz']);
% save_untouch_nii(v2,[velDir '/' velName '-bkCorr-2.nii.gz']);

save_untouch_nii(v0,[velDir '/' velName '-polyCorr-0.nii.gz']);
save_untouch_nii(v1,[velDir '/' velName '-polyCorr-1.nii.gz']);
save_untouch_nii(v2,[velDir '/' velName '-polyCorr-2.nii.gz']);

disp('Saved polynomial corrected velocity volumes ...');

%% Save blood pool mask
cd(reconDir);

bp_nii = v0;
bp_nii.img = logical(BP);
save_untouch_nii(bp_nii,'mask/mask_blood_pool.nii.gz');


% fn end
end