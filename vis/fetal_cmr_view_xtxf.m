
function [xtRcn, xfRcn] = fetal_cmr_view_xtxf( xtRcn, sliceNum, xfLine, cLims )

% Functions
dimT = 3;
xf2xt = @( xf ) ifft( ifftshift( xf, dimT ), [], dimT );
xt2xf = @( xt ) fftshift( fft( xt, [], dimT ), dimT );

xtRcn = squeeze( xtRcn );
xfRcn = xt2xf( xtRcn );

 
implay_RR(abs(xtRcn(:,:,:,sliceNum)));

figure;
set(gcf,'Position', [900   300   1000   700]);
subplot(1,2,1);
imagesc(abs(xtRcn(:,:,1,sliceNum))); 
% axis('image'); 
set(gca,'YLimMode','auto'); 
colormap('gray');
title( ['x-t - Slice No: ' num2str(sliceNum)] );
freezeColors;

if nargin == 3
    
    cLims(1) = 0;
    cLims(2) = 5 * mean(abs(xfRcn(:)));

    subplot(1,2,2);
    imagesc( squeeze(abs(xfRcn(:,xfLine,:,sliceNum))), [cLims(1), cLims(2)] );
    colormap('parula'); colorbar;
%     axis('image');
    set(gca,'YLimMode','auto'); 
    title( ['x-f - Slice No: ' num2str(sliceNum)] );
    
elseif nargin > 3    

    subplot(1,2,2);
    imagesc( squeeze(abs(xfRcn(:,xfLine,:,sliceNum))), [cLims(1), cLims(2)] );
    colormap('parula'); colorbar;
%     axis('image');
    set(gca,'YLimMode','auto'); 
    title( ['x-f - Slice No: ' num2str(sliceNum)] );

end









