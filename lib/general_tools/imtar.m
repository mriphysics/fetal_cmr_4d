function imtar(Image,Cl,Cu,cmapcolour)

% imtar
%
% Quick function so I don't keep having to type out a long boring line!
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     disp(' I ? Baz');
    if nargin == 1
        figure; imagesc(Image(:,:)), colormap(gray), colorbar, axis image;
    elseif nargin > 1 && nargin < 4
        figure, imagesc(Image(:,:),[Cl Cu]), colormap(gray), colorbar, axis image;
    elseif nargin == 4
        figure, imagesc(Image(:,:),[Cl Cu]), colormap(cmapcolour), colorbar, axis image;
    end