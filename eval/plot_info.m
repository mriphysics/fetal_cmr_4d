function [I,H] = plot_info( infoFile )
% PLOT_INFO  read and plot data from info.tsv file output from reconstructionCardiac
% 
%   [I,H] = PLOT_INFO( infoFile ) reads a .tsv file, plots the data and 
%   returns structured data in I and a cell array of figure handles in H

%   jfpva (joshua.vanamerom@kcl.ac.uk) 


%% Get Data

I = read_info_tsv( infoFile );


%% Pre-Process

stackID = unique( I.StackIndex );

stackDesc = cell( size( stackID ) );

sliceID   = cell( size( stackID ) );

for iSt = 1:numel(stackID)
    
    % stack description
    stackFilePath = I.File(find(I.StackIndex==stackID(iSt),1));  
    [ ~, desc, ~ ] = fileparts( stackFilePath.char );
    while any( desc == '.' )
        [ ~, desc, ~ ] = fileparts( desc );
    end
    stackDesc{iSt} = desc;
    
    % slice-location IDs
    sliceID{iSt} = unique( I.StackLocIndex( I.StackIndex == stackID(iSt) ) );

end

% Figures 
SubFigSizeX = 150;
SubFigSizeY = 125;

setlocation = @( h ) setappdata( h, 'SubplotDefaultAxesLocation', [0.05, 0.025, 0.925, 0.95] );
setinset    = @( h ) setappdata( h, 'SubplotDefaultInset', [0, 0, 0, 0]);

H = cell( size( stackID ) );


%% Plot

nRow = 5;

limTranslation  = 10;
limRotation     = 10;
limDisplacement = 10;
limTRE          = limDisplacement;

cLines          = lines(12);
cWeight         = cLines(5,:);
cDisp           = cLines(6,:);
cDispWeighted   = cLines(12,:);
cTRE            = cLines(4,:);

for iSt = 1:numel(stackID)

    nSlice = numel(sliceID{iSt}) + 1;
    
    hFig = figure('Name',[stackDesc{iSt},'_transformations_v_time'],'Position',[0,250,nSlice*SubFigSizeX,nRow*SubFigSizeY]);
    
    setlocation( hFig );
    setinset( hFig );
    
    for iSl = 1:(nSlice-1)
        ind = I.StackIndex == stackID(iSt) & I.StackLocIndex == sliceID{iSt}(iSl);
        t   = I.Time( ind );
        w   = I.Weight( ind );
        tx  = I.TranslationX( ind );
        ty  = I.TranslationY( ind );
        tz  = I.TranslationZ( ind );
        rx  = I.RotationX( ind );
        ry  = I.RotationY( ind );
        rz  = I.RotationZ( ind );
        dx  = I.MeanDisplacementX( ind );
        dy  = I.MeanDisplacementY( ind );
        dz  = I.MeanDisplacementZ( ind );
        d   = I.MeanDisplacement( ind );
        dw  = I.WeightedMeanDisplacement( ind );
        tre = I.TRE( ind );
        iRow = 1;
        % Translation
        subplot(nRow,nSlice,(iRow-1)*nSlice+iSl) 
        plot(t,tx,t,ty,t,tz,'LineWidth',1)
        hAx = gca;
        hAx.XLim = [min(t),max(t)];
        hAx.XTick = 0:2:max(t);
        hAx.XTickLabel = '';
        hAx.YLim = limTranslation*[-1,+1];
        hAx.YTick = hAx.YLim(2)*[-1,-0.5,0,0.5,1];
        grid on
        if iRow == 1 
            title( sprintf( 'slice %i', sliceID{iSt}(iSl)+1 ) )
        end
        if iSl == 1
            ylabel( 'Translation (mm)' )
        else
            hAx.YTickLabel = '';
        end
        if iSl == nSlice-1
            hLgnd = legend( {'tx','ty','tz'}, 'Orientation', 'vertical', 'Location', 'West' );
            hAx = subplot(nRow,nSlice,(iRow-1)*nSlice+iSl+1);
            axis off
            hLgnd.Position(1) = hAx.Position(1);
        end
        iRow = iRow + 1;
        % Rotation
        subplot(nRow,nSlice,(iRow-1)*nSlice+iSl)
        plot(t,rx,t,ry,t,rz,'LineWidth',1)
        hAx = gca;
        hAx.XLim = [min(t),max(t)];
        hAx.XTick = 0:2:max(t);
        hAx.XTickLabel = '';
        hAx.YLim = limRotation*[-1,+1];
        hAx.YTick = hAx.YLim(2)*[-1,-0.5,0,0.5,1];
        grid on
        if iRow == 1 
            title( sprintf( 'slice %i', sliceID{iSt}(iSl)+1 ) )
        end
        if iSl == 1
            ylabel( 'Rotation (deg.)' )
        else
            hAx.YTickLabel = '';
        end
        if iSl == nSlice-1
            hLgnd = legend( {'rx','ry','rz'}, 'Orientation', 'vertical', 'Location', 'West' );
            hAx = subplot(nRow,nSlice,(iRow-1)*nSlice+iSl+1);
            axis off
            hLgnd.Position(1) = hAx.Position(1);
        end
        iRow = iRow + 1;
        % Displacement
        subplot(nRow,nSlice,(iRow-1)*nSlice+iSl) 
        plot(t,dx,t,dy,t,dz,'LineWidth',1)
        hold on
        hAx = gca;
        hAx.XLim = [min(t),max(t)];
        hAx.XTick = 0:2:max(t);
        hAx.XTickLabel = '';
        hAx.YLim = limDisplacement*[-1,+1];
        hAx.YTick = hAx.YLim(2)*[-1,-0.5,0,0.5,1];
        grid on
        if iRow == 1 
            title( sprintf( 'slice %i', sliceID{iSt}(iSl)+1 ) )
        end
        if iSl == 1
            ylabel( 'Displacement (mm)' )
        else
            hAx.YTickLabel = '';
        end
        if iSl == nSlice-1
            hLgnd = legend( {'dx','dy','dz'}, 'Orientation', 'vertical', 'Location', 'West' );
            hAx = subplot(nRow,nSlice,(iRow-1)*nSlice+iSl+1);
            axis off
            hLgnd.Position(1) = hAx.Position(1);
        end
        iRow = iRow + 1;
        % TRE
        subplot(nRow,nSlice,(iRow-1)*nSlice+iSl) 
        hold on
        plot(t,d,'Color',cDisp,'LineWidth',1)
        plot(t,dw,'--','Color',cDispWeighted,'LineWidth',0.5)
        plot(t,tre,'Color',cTRE,'LineWidth',1)
        hold off
        hAx = gca;
        hAx.XLim = [min(t),max(t)];
        hAx.XTick = 0:2:max(t);
        hAx.XTickLabel = '';
        hAx.YLim = limTRE*[0,1];
        hAx.YTick = hAx.YLim(2)*(0:0.5:1);
        grid on
        box on
        if iRow == 1 
            title( sprintf( 'slice %i', sliceID{iSt}(iSl)+1 ) )
        end
        if iSl == 1
            ylabel( 'Disp./TRE (mm)' )
        else
            hAx.YTickLabel = '';
        end
        if iSl == nSlice-1
            hLgnd = legend( { 'd','d_{weighted}','TRE' }, 'Orientation', 'vertical', 'Location', 'West' );
            hAx = subplot(nRow,nSlice,(iRow-1)*nSlice+iSl+1);
            axis off
            hLgnd.Position(1) = hAx.Position(1);
        end
        iRow = iRow + 1;
        % Weight
        subplot(nRow,nSlice,(iRow-1)*nSlice+iSl) 
        plot(t,w,'Color',cWeight,'LineWidth',1.5)
        hAx = gca;
        hAx.XLim = [min(t),max(t)];
        hAx.XTick = 0:2:max(t);
        hAx.XTickLabel = '';
        hAx.YLim = [0,1];
        hAx.YTick = hAx.YLim(2)*(0:0.25:1);
        hAx.YTickLabel = {'0','','0.5','0','1'};
        grid on
        if iRow == 1 
            title( sprintf( 'slice %i', sliceID{iSt}(iSl)+1 ) )
        end
        if iSl == 1
            ylabel( 'p_k^{frame}' )
        else
            hAx.YTickLabel = '';
        end
        if iSl == nSlice-1
            hLgnd = legend( {'p_k^{frame}'}, 'Orientation', 'vertical', 'Location', 'West' );
            hAx = subplot(nRow,nSlice,(iRow-1)*nSlice+iSl+1);
            axis off
            hLgnd.Position(1) = hAx.Position(1);
        end
        iRow = iRow + 1;
    end
    
    H{iSt} = hFig;
    
end


end  % plot_info(...)
