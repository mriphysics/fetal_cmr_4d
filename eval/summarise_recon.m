function I = summarise_recon( recondir, cardsyncdir, varargin )
%SUMMARISE_RECON  summarises volumetric cine reconstruction results
%
%   I = SUMMARISE_RECON( recondir, cardsyncdir ) returns results in
%   structure I given directories for cine reconstruction and cardiac
%   synchronisation data.
% 
%   SUMMARISE_RECON( ..., 'name', value ) specifies optional input 
%   arguments. See code for name-value pairs.

%   jfpva (joshua.vanamerom@kcl.ac.uk) 


%% Optional Input Argument Default Values

default.infoFileName     = 'info.tsv';
default.logMainFileName  = 'log-main.txt';
default.logReconFileName = 'log-reconstruction.txt';
default.cardsyncFileName = 'results_interslice_cardsync.mat';
default.niftiPath        = fullfile(fileparts(mfilename('fullpath')),'..','lib','nifti');  % path to NIfTI toolbox
default.irtkPath         = '/usr/local/build/bin/';  % path to IRTK bin 
default.isVerbose        = false;


%% Parse Input

p = inputParser;

if  verLessThan('matlab','8.2')
    add_param_fn = @( parseobj, argname, defaultval, validator ) addParamValue( parseobj, argname, defaultval, validator );
else
    add_param_fn = @( parseobj, argname, defaultval, validator ) addParameter( parseobj, argname, defaultval, validator );
end

addRequired(   p, 'recondir', ... 
    @(x) validateattributes( x, {'char'}, {'vector'}, mfilename) );

addRequired(   p, 'cardsyncdir', ... 
    @(x) validateattributes( x, {'char'}, {'vector'}, mfilename) );

add_param_fn(   p, 'infofilename', default.infoFileName, ...
    @(x) validateattributes( x, {'char'}, {'vector'}, mfilename) );

add_param_fn(   p, 'logmainfilename', default.logMainFileName, ...
    @(x) validateattributes( x, {'char'}, {'vector'}, mfilename) );

add_param_fn(   p, 'logreconfilename', default.logReconFileName, ...
    @(x) validateattributes( x, {'char'}, {'vector'}, mfilename) );

add_param_fn(   p, 'cardsyncfilename', default.cardsyncFileName, ...
    @(x) validateattributes( x, {'char'}, {'vector'}, mfilename) );

add_param_fn(   p, 'niftipath', default.niftiPath, ...
    @(x) validateattributes( x, {'char'}, {'vector'}, mfilename) );

add_param_fn(   p, 'irtkpath', default.irtkPath, ...
    @(x) validateattributes( x, {'char'}, {'vector'}, mfilename) );

add_param_fn(   p, 'verbose', default.isVerbose, ...
    @(x) validateattributes( x, {'logical'}, {'scalar'}, mfilename) );

parse( p, recondir, cardsyncdir, varargin{:} );

reconDir            = p.Results.recondir;
cardsyncDir         = p.Results.cardsyncdir;
infoFileName        = p.Results.infofilename;
logMainFileName     = p.Results.logmainfilename;
logReconFileName    = p.Results.logreconfilename;
cardsyncFileName    = p.Results.cardsyncfilename;
niftiPath           = p.Results.niftipath;
irtkPath            = p.Results.irtkpath;
isVerbose           = p.Results.verbose;


%% Dependencies

% Temporarily add Matlab toolboxes
origPath  = path;
resetPath = onCleanup( @() path( origPath ) );
addpath( genpath( niftiPath ) )

% Temporarily add IRTK to matlab environment
origEnv  = getenv('PATH');
resetEnv = onCleanup( @() setenv( 'PATH', origEnv ) );
setenv( 'PATH', [ origEnv ':' irtkPath ] );  


%% File Paths

infoFilePath    = fullfile( reconDir, infoFileName );
logMainFilePath = fullfile( reconDir, logMainFileName );
logReconFilePath = fullfile( reconDir, logReconFileName );
cardsyncFilePath = fullfile( cardsyncDir, cardsyncFileName );


%% Load Info 

T = read_info_tsv( infoFilePath, logReconFilePath );


%% Load Cardsync Data

M = matfile( cardsyncFilePath );
R = M.I;


%% Init Output

I = struct;
S = struct;


%% Identify Outside Slices

outsideImages           = T.InputIndex( T.Outside==1 ) + 1;
smallImages             = T.InputIndex( T.Small==1 ) + 1;
indNoOverlap            = union( outsideImages, smallImages );  % indices considered outside
indOverlap              = setdiff( T.InputIndex + 1, indNoOverlap );
outsideLocIndex         = T.LocIndex( indNoOverlap );
outsideStackDynIndex    = T.StackDynIndex( indNoOverlap );
outsideSlices           = [];
for locIndex = unique( outsideLocIndex )'
    nFrameDyn = max(T.StackDynIndex(T.LocIndex==locIndex))+1;
    if isequal( (1:nFrameDyn)-1, sort( outsideStackDynIndex( outsideLocIndex == locIndex ) )' )
        outsideSlices = cat( 2, outsideSlices, locIndex + 1 );
    end
end
insideSlices            = setdiff( T.LocIndex+1, outsideSlices );
S.inside = false( 1, max(T.LocIndex)+1 );
S.inside( insideSlices) = true;

T.Inside = ~( T.Outside | T.Small );


%% Calculate Heart Rate Results

tRRused = cell2mat( [ R.tRR ] );
tRRest  = cell2mat( [ R.tRRest ] );
tRRest( outsideSlices ) = NaN;
tRR     = tRRused;
tRR( tRRest ~= tRRused ) = NaN;
S.heartrate.tRR   = tRR;
S.heartrate.mean  = mean( tRR, 'omitnan' );
S.heartrate.std   = std( tRR, 'omitnan' );


%% Determine Number of Reconstruction Iterations

DE = dir( fullfile( reconDir, 'sr_iterations', 'errorstack000_mc*sr*.nii.gz' ) );

nmc = nan(1,numel(DE));
nsr = nan(1,numel(DE));

for iF = 1:numel(DE)
    nmc(iF) = str2double(DE(iF).name(17:18));
    nsr(iF) = str2double(DE(iF).name(21:22));
end

Nmc = max(nmc);
Nsr = max(nsr);


%% Calculate Outlier Results v2

outlier.numincludedim     = sum( T.Included & T.Inside );
outlier.numexcludedim     = sum( T.Excluded & T.Inside );
outlier.percentexcludedim = 100 * outlier.numexcludedim / ( outlier.numincludedim + outlier.numexcludedim );

DW = dir( fullfile( reconDir, 'weight*.nii.gz' ) );
DV = dir( fullfile( reconDir, 'sr_iterations', sprintf( 'simweightstack*_mc%02isr%02i.nii.gz', Nmc, 0 ) ) );
DE = dir( fullfile( reconDir, 'sr_iterations', sprintf( 'errorstack*_mc%02isr%02i.nii.gz', Nmc, Nsr ) ) );
DM = dir( fullfile( reconDir, 'stack*.nii.gz' ) );

pVoxel = [];
eVoxel = [];
outlier.pvox = cell( 1, max(T.InputIndex)+1 );
eVoxIm = cell( 1, max(T.StackIndex)+1 );
pVoxIm = cell( 1, max(T.StackIndex)+1 );
maskIm = cell( 1, max(T.StackIndex)+1 );
for iF = 1:numel(DW)
    W = load_untouch_nii( fullfile( DW(iF).folder, DW(iF).name ) );
    E = load_untouch_nii( fullfile( DE(iF).folder, DE(iF).name ) );
    V = load_untouch_nii( fullfile( DV(iF).folder, DV(iF).name ) );
    M = load_untouch_nii( fullfile( DM(iF).folder, DM(iF).name ) );
    eVoxIm{iF} = E.img;
    pVoxIm{iF} = W.img;
    maskIm{iF} = ( (M.img~=-1) & (V.img>0.99) );
    eVoxel = [ eVoxel; eVoxIm{iF}(maskIm{iF}) ];
end
T.NumVoxel = zeros( size( T.InputIndex ) );
for i = T.InputIndex(:)'+1
    indStack = T.StackIndex(i)+1;
    indSlice = T.StackLocIndex(i)+1;
    indDyn   = T.StackDynIndex(i)+1;
    p = pVoxIm{indStack}(:,:,indSlice,indDyn);
    m = maskIm{indStack}(:,:,indSlice,indDyn);
    outlier.pvox{i} = p(m);
    T.NumVoxel(i) = sum(m(:));
    if any( i == indOverlap )
        pVoxel = [ pVoxel; outlier.pvox{i} ];
    end
end

outlier.evoxinside = eVoxel;
outlier.pvoxinside = pVoxel;
outlier.percentexcludedvox = 100 * sum(pVoxel<0.5) / numel(pVoxel);

logFilename = fullfile( reconDir, 'log-reconstruction.txt' );
if exist( 'logFilename', 'var' )
    cmd = sprintf( 'grep "Voxel-wise robust statistics parameters:" %s | tail -1', logFilename );
    [ ~, result ] = system( cmd );
    if ~isempty( result )
        C = strsplit(result,' ');
        outlier.voxel.inclass.std       = str2double( C{7} );
        outlier.voxel.inclass.mix       = str2double( C{10} );
        outlier.voxel.outclass.density  = str2double( C{13} );
    else
        outlier.voxel.inclass.std       = NaN;
        outlier.voxel.inclass.mix       = NaN;
        outlier.voxel.outclass.density  = NaN;
    end
end


logFilename = fullfile( reconDir, 'log-reconstruction.txt' );
if exist( 'logFilename', 'var' )
    cmd = sprintf( 'grep "Slice robust statistics parameters:" %s | tail -1', logFilename );
    [ ~, result ] = system( cmd );
    if ~isempty( result )
        C = strsplit(result,' ');
        outlier.frame.inclass.mean   = str2double( C{6} );
        outlier.frame.outclass.mean  = str2double( C{7} );
        outlier.frame.inclass.std    = str2double( C{9} );
        outlier.frame.outclass.std   = str2double( C{10} );
        outlier.frame.inclass.mix    = str2double( C{12} );
        outlier.frame.outclass.mix   = str2double( C{13} );
    else
        outlier.frame.inclass.mean   = NaN;
        outlier.frame.outclass.mean  = NaN;
        outlier.frame.inclass.std    = NaN;
        outlier.frame.outclass.std   = NaN;
        outlier.frame.inclass.mix    = NaN;
        outlier.frame.outclass.mix   = NaN;
    end
end

I.outlier = outlier;


%% Plot Outlier Rejection

minPotential = 0.01;

rgbIn   = [127,255,79]/255;
rgbOut = [215,25,28]/255;  
rgbProb = [253,174,97]/255; 
rgbMix  = [117,158,79]/255;

tmp = bone(5);
rgbBar = tmp(4,:);

if (isVerbose) 

    hFig = figure( 'Name', 'outlier_rejection' );
    hFig.Position(4) = hFig.Position(4)/2;
        
    % Voxel-Wise 
    subplot(1,2,1)
    yyaxis left
    xLimVox = prctile(abs(outlier.evoxinside),99.8)*[-1,+1];
    xVox = linspace(xLimVox(1),xLimVox(2),1000);
    if isnan( outlier.voxel.inclass.mix )
        inVoxelPdf  = nan(size(xVox));
        outVoxelPdf = nan(size(xVox));
        voxProb     = nan(size(xVox));
    else
        inVoxelPdf  = outlier.voxel.inclass.mix *  normpdf( xVox, 0,  outlier.voxel.inclass.std );
        outVoxelPdf = ( 1 - outlier.voxel.inclass.mix ) * outlier.voxel.outclass.density * ones( size( xVox ) );
        voxProb = inVoxelPdf ./ ( inVoxelPdf + outVoxelPdf );
    end
    binWidthVox = 0.0373 * range(xLimVox); %28;
    binEdgesVox = [-flip((binWidthVox/2):binWidthVox:(-xLimVox(1)+binWidthVox/2)),(binWidthVox/2):binWidthVox:(xLimVox(2)+binWidthVox/2)];
    hHistVox = histogram(outlier.evoxinside,binEdgesVox,'Normalization','pdf','FaceColor',rgbBar);
    hold on
    plot(xVox,inVoxelPdf+outVoxelPdf,'-','LineWidth',3,'Color',rgbMix);
    plot(xVox,outVoxelPdf,'--','LineWidth',1.5,'Color',rgbOut);
    plot(xVox,inVoxelPdf,'--','LineWidth',1.5,'Color',rgbIn);
    yyaxis right
    plot(xVox,voxProb,'LineWidth',2.5,'Color',rgbProb)
    yyaxis left
    yLimVox = [0,max(max(hHistVox.Values),max(inVoxelPdf+outVoxelPdf))];
    set(gca,'XLim',xLimVox,'YLim',yLimVox)
    set(gca,'YColor',[0,0,0],'YTickLabel',{})
    ylabel('P(\it{e})','FontSize',14)
    xlabel('\it{e}','FontSize',14)
    yyaxis right
    set(gca,'YColor',[0,0,0],'YTickLabel',{})
    ylabel('\it{p^{voxel}}','FontSize',14)    
        
    % Frame-Wise 
    subplot(1,2,2) 
    xLimFrm = [0,min(1,prctile(T.Potential(T.Potential>=minPotential),99.8))];
    xFrm = linspace(xLimFrm(1),xLimFrm(2),1000);
    if isnan( outlier.voxel.inclass.mix )
        inFramePdf  = nan(size(xFrm));
        outFramePdf = nan(size(xFrm));
        frameProb   = nan(size(xFrm));
    else
        inFramePdf  = outlier.frame.inclass.mix *  normpdf( xFrm, outlier.frame.inclass.mean,  outlier.frame.inclass.std );
        outFramePdf = outlier.frame.outclass.mix * normpdf( xFrm, outlier.frame.outclass.mean, outlier.frame.outclass.std );
        frameProb   = inFramePdf ./ ( inFramePdf + outFramePdf );
    end
    binWidthFrm = 0.0571 * range(xLimFrm); %0.02;
    binEdgesFrm = xLimFrm(1):binWidthFrm:(xLimFrm(2)+binWidthFrm/2);
    yyaxis left
    hHistFrm = histogram(T.Potential(T.Potential>=minPotential),binEdgesFrm,'Normalization','pdf','FaceColor',rgbBar);
    hold on
    plot(xFrm,inFramePdf+outFramePdf,'-','LineWidth',3,'Color',rgbMix);
    plot(xFrm,outFramePdf,'--','LineWidth',1.5,'Color',rgbOut);
    plot(xFrm,inFramePdf,'--','LineWidth',1.5,'Color',rgbIn);
    yyaxis right
    plot(xFrm,frameProb,'LineWidth',2.5,'Color',rgbProb)
    yyaxis left
    yLimFrm = [0,max(max(max(hHistFrm.Values),1.05*max(inFramePdf+outFramePdf)),eps)];
    set(gca,'XLim',xLimFrm,'YLim',yLimFrm)
    set(gca,'YColor',[0,0,0],'YTickLabel',{})
    ylabel('P(\it{q})','FontSize',14)
    xlabel('\it{q}','FontSize',14)
    yyaxis right
    set(gca,'YColor',[0,0,0],'YLim',[0,1],'YTick',[0,1],'YTickLabel',{})
    ylabel('\it{p}^{frame}','FontSize',14)
    
    % NOTE: issue following from commandline to get liklihood symbol matching notation in manuscript
    %{
          subplot(1,2,1);
          hAx1 = gca;
          hAx1.YAxis(1).Label.String = '?(\it{e})';
          subplot(1,2,2);
          hAx2 = gca;
          hAx2.YAxis(1).Label.String = '?(\it{q})';
    %}
    
end


%% Summarise Transformations

I.transformations = summarise_transformations( reconDir, T, 'niftipath', niftiPath, 'irtkpath', irtkPath );


%% Get Displacement Results Calculated During Cine Volume Reconstruction

cmd = sprintf( 'grep "Mean Displacement:" %s | awk -F''[ ]'' ''{print $6}''', logMainFilePath );
[ ~, result ] = system( cmd );
I.disp.allvox.mean = str2double( result );
I.disp.insidevox.stdmean = std( T.MeanDisplacement(T.Inside), 'omitnan' );
I.disp.insidevox.mean = sum( T.NumVoxel(T.Inside) .* T.MeanDisplacement(T.Inside) ) / sum( T.NumVoxel(T.Inside) ); 

S.disp.mean = nan( 1, max( T.LocIndex )+1 );
S.disp.std  = nan( 1, max( T.LocIndex )+1 );
locIndices = unique( T.LocIndex )'; 
for locIndex = locIndices
    S.disp.mean(locIndex+1) = mean(T.MeanDisplacement(T.LocIndex==locIndex),'omitnan'); 
    S.disp.std(locIndex+1)  = std( T.MeanDisplacement(T.LocIndex==locIndex),'omitnan'); 
end
S.disp.insideim.meanstd = mean( S.disp.std(insideSlices), 'omitnan' );


%% Combine Results

I.slice = S;
I.T = T;


end  % summarise_recon(...)