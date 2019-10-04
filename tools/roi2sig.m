function [Sig, meanSig, ROI] = roi2sig(im,roipoly_slice,ROI)


% im - must be 3D. 3rd dimension is b-shells etc.

if nargin == 1
    roipoly_slice = 1;
    ROI = [];
end

% if nargin == 2
%     ROI = [];
% end


%% Draw ROI
clear Sig

if nargin == 3
    % ROI supplied to function
else
%     imtar(im(:,:,roipoly_slice));
%     ROI = double(roipoly); close;
    ROI = freeROI(im(:,:,roipoly_slice));
    close
end

Sig = zeros(length(find(ROI)),size(im,3));

for i = 1:size(im,3)   
    tempIM = ROI .* double(im(:,:,i));
    Sig(:,i) = tempIM(find(ROI));
end

meanSig = mean(Sig,1);

figure;
plot(meanSig,'o');


end
%fn end