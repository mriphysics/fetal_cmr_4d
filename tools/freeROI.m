function [ mask ] = freeROI( im, active, iterations )
% Laurence Jackson, CABI, UCL, 2015
% 
% Draw freehand region of interest on image im
% im            = image to draw on
% mask          = mask of region
% active        = 0 - (defualt) no active contours
%                 1 - Uses active contouring to find edges
% iterations    = number of iterations for active contouring (default = 10)
%%
figure('units','normalized','outerposition',[0 0 1 1])
imagesc(im); colormap('gray'); axis image;
title('Draw freehand:');
h = imfreehand(gca);
mask = createMask(h);

    imagesc(im); colormap('gray'); axis image; 
    hold on;
    [b] = bwboundaries(mask);
    b = b{:};
    plot(b(:,2), b(:,1), 'w', 'LineWidth', 2)
    pause(1);
    
% close h;
%% active contour
if nargin > 1 && active == 1
    if nargin > 2 
        it = iterations;
    else
        it = 10;        
    end
    
    mask2 = activecontour(im, mask, it,  'edge' );
    
    imagesc(im); hold on;
    [b] = bwboundaries(mask2);
    b = b{:};
    plot(b(:,2), b(:,1), 'w', 'LineWidth', 2)  
    mask = mask2;
    
end

return 


