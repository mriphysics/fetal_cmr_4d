# fetal_cmr_4d

fetal whole-heart 4d reconstruction using motion-corrected multi-planar real-time MRI

## Publications

__Fetal whole-heart 4D imaging using motion-corrected multi-planar real-time MRI.__  
Joshua FP van Amerom, David FA Lloyd, Maria Deprez, Anthony N Price, Shaihan J Malik, Kuberan Pushparajah, Milou PM van Poppel, Mary A Rutherford, Reza Razavi, Joseph V Hajnal.  
13 Apr 2019. Magnetic Resonance in Medicine. 2019. doi: [10.1002/mrm.27798](https://doi.org/10.1002/mrm.27798) (_accepted, peer-reviewed+revised_)  
05 Dec 2018. arXiv: [1812.02249](https://arxiv.org/abs/1812.02249) (_preprint_)  

## Directories

__4drecon__ - preprocessing and 4D reconstruction scripts  

__cardsync__ - cardiac synchronisation   

__eval__ - summarise and evaluate results

__irtk_cardiac4d__ - 4D reconstruction submodule linked to [github.com/jfpva/irtk_cardiac4d](https://github.com/jfpva/irtk_cardiac4d), built on the Image Registration Toolkit (IRTK) v2.0.  

__ktrecon__ - k-t sense reconstruction submodule linked to  [github.com/jfpva/ktrecon](https://github.com/jfpva/ktrecon)  


## Installation

Build instructions for irtk can be found at [sites.google.com/site/mariakuklisovamurgasova/software](https://sites.google.com/site/mariakuklisovamurgasova/software).


## External Dependencies

__ReconFrame__ - software platform providing the tools and the functionality to develop and execute a complete image reconstruction of Philips MR data ([gyrotools.com/gt/index.php/products/reconframe](https://www.gyrotools.com/gt/index.php/products/reconframe))  


## Framework 

Framework for 4D cine reconstruction, consists of:

1. acquisition and reconstruction of multi-planar dynamic2DMRI; 
2. an initial motion correction stage to achieve rough spatial alignment of the fetal heart using temporal mean (i.e., static) images for stack-stack registration followed by slice-volume registration interleaved with static volume (3D) reconstruction; 
3. cardiac synchronisation, including heart rate estimation and slice-slice cardiac cycle alignment; and 
4. further motion-correction using dynamic image frames interleaved with 4D reconstruction; and 
5. a final 4D cine reconstruction, including outlier rejection.  

![](4d_framework.png)  

User-specified ROIs, and identification of a target stack for stack-stack registration are the only manual preparations required for reconstruction.


## Steps

__Setup__  \
create working directories,  \
e.g., in shell:  
```shell 
RECONDIR=~/path/to/recon/directory
mkdir $RECONDIR
mkdir $RECONDIR/data
mkdir $RECONDIR/ktrecon
mkdir $RECONDIR/mask
mkdir $RECONDIR/cardsync
```

1. __MRI__
    - acquire 2D multi-planar real-time MRI data
    - reconstruct images using `ktrecon`, \
    e.g., for each stack, in Matlab:  
        ```matlab
        reconDirPath        = '~/path/to/recon/directory';
        seriesNo            = 0;
        rawDataFilePath     = '~/path/to/rawdata.lab';
        senseRefFilePath    = '~/path/to/senserefscan.lab';
        coilSurveyFilePath  = '~/path/to/coilsurveyscan.lab';
        outputDirPath       = fullfile( reconDirPath, 'ktrecon' );
        outputStr           = sprintf( 's%02i', seriesNo );
        reconOpts           = { 'GeometryCorrection', 'Yes' };

        mrecon_kt(   rawDataFilePath, ...
                    'senseref', senseRefFilePath, ...
                    'coilsurvey', coilSurveyFilePath, ...
                    'outputdir', outputDirPath, ...
                    'outputname', outputStr, ...
                    'patchversion', patchVersion,...
                    'reconoptionpairs', reconOpts )
        ```
    - further processsing
        - copy/move all magnitude-valued DC (s\*_dc_ab.nii.gz) and real-time (s\*_rlt_ab.nii.gz) files from 'ktrecon' directory to 'data' directory \
        e.g., in shell: 
            ```shell
            cp ktrecon/s*_dc_ab.nii.gz data;
            cp ktrecon/s*_rlt_ab.nii.gz data;
            ```
        - manually draw fetal heart masks for each `sXX_dc_ab.nii.gz` file (e.g., using the [Medical Imaging ToolKit (MITK) Workbench](http://mitk.org/wiki/Downloads#MITK_Workbench))
            - draw ROI containing fetal heart and great vessels for each slice
            - save segmentation as `sXX_mask_heart.nii.gz` segmentation in 'mask' directory
        - run `preproc` in Matlab,
            ```matlab
            reconDir = '~/path/to/recon/directory';
            S = preproc( reconDir );
            save( fullfile( reconDir, 'data', 'results.mat' ), 'S', '-v7.3' );
            ```
        - optionally, manually specify
            - target stack by changing value in 'data/tgt_stack_no.txt' (stacks are index 1,2,...)
            - excluded stacks/slices/frames by specifying in 'data/force_exclude_*.txt' (stacks/slices/frames are zero-indexed)
2. __Motion-Correction (static)__
    - create 3D mask of fetal chest
        - recon reference volume, \
        e.g., in shell: 
        ```shell
        RECONDIR=~/path/to/recon/directory
        ./recon_ref_vol.bash $RECONDIR ref_vol
        ```
        - draw fetal chest ROI using 'ref_vol.nii.gz' as a reference (e.g., using MITK)
        - save segmentation to 'mask' directory as  'mask_chest.nii.gz'
            - _note:_ the orientation of all later 3D/4D reconstructions is determined by this mask file; the orientation can be changed by applying a transformation to 'mask_chest.nii.gz' prior to further reconstructions
    - static (slice-wise) motion-correction, \
    e.g., in shell: 
        ```shell
        RECONDIR=~/path/to/recon/directory˜
        ./recon_dc_vol.bash $RECONDIR dc_vol
        ```
3. __Cardiac Synchronisation__
    - heart-rate estimation
        - run `cardsync_intraslice`, in Matlab:
            ```matlab
            reconDir    = '~/path/to/recon/directory';
            dataDir     = fullfile( reconDir, 'data' );
            cardsyncDir = fullfile( reconDir, 'cardsync' );
            M = matfile( fullfile( dataDir, 'results.mat' ) );
            S = cardsync_intraslice( M.S, 'resultsDir', cardsyncDir, 'verbose', true );
            ```
    - slice-slice synchronisation
        - recon cine volume for each slice, \
        e.g., in shell: 
            ```shell
            RECONDIR=~/path/to/recon/directory
            ./recon_slice_cine.bash $RECONDIR
            ```
        - optionally, specify target slice by creating file 'data/tgt_slice_no.txt' containing target slice number (indexed starting at 1)
        - run `cardsync_interslice`, in Matlab:
            ```matlab
            % setup
            reconDir    = '~/path/to/recon/directory';
            dataDir     = fullfile( reconDir, 'data' );
            cardsyncDir = fullfile( reconDir, 'cardsync' );
            cineDir     = fullfile( reconDir, 'slice_cine_vol' );    
            M = matfile( fullfile( cardsyncDir, 'results_cardsync_intraslice.mat' ) );
            
            % target slice
            tgtLoc = NaN;
            tgtLocFile = fullfile( dataDir, 'tgt_slice_no.txt' );
            if exist( tgtLocFile , 'file' )
              fid = fopen( tgtLocFile, 'r' );
              tgtLoc = fscanf( fid, '%f' );
              fclose( fid );
            end
            
            % excluded slices
            excludeSlice = [];
            excludeSliceFile = fullfile( dataDir, 'force_exclude_slice.txt' );
            if exist( excludeSliceFile , 'file' )
              fid = fopen( excludeSliceFile, 'r' );
              excludeSlice = fscanf( fid, '%f' ) + 1;  % NOTE: slice locations in input file are zero-indexed
              fclose( fid );
            end
            
            % slice-slice cardiac synchronisation
            S = cardsync_interslice( M.S, 'recondir', cineDir, 'resultsdir', cardsyncDir, 'tgtloc', tgtLoc, 'excludeloc', excludeSlice );
            ```
4. __Motion-Correction (dynamic)__
    - performed interleaved with 4D Reconstruction
5. __4D Volumetric Reconstruction__
    - recon 4D (cine) volume, \
    e.g., in shell: 
    ```shell
    RECONDIR=~/path/to/recon/directory
    ./recon_cine_vol.bash $RECONDIR cine_vol
    ```

__Summarise__  \
e.g., in Matlab:  
```matlab
S = summarise_recon( '~/path/to/recon/directory/cine_vol', '~/path/to/recon/directory/cardsync', 'verbose', true );
I = plot_info( '~/path/to/recon/directory/cine_vol/info.tsv');
```