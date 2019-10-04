function IMrot = rot_ms(IM,deg)

%% rot_ms - rotate multi-slice data
%
% Quick function to extend rot90 to 3- or 4-dimensional data.
%
% - Tom Roberts
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% IMrot = zeros(size(IM));

if numel(size(IM)) == 3

    dim3 = size(IM,3);
    
    for i = 1:dim3;
            IMrot(:,:,i) = rot90(IM(:,:,i),deg);
    end
    
elseif numel(size(IM)) == 4

    dim3 = size(IM,3);
    dim4 = size(IM,4);

    for i = 1:dim3;
        for j = 1:dim4;
            IMrot(:,:,i,j) = rot90(IM(:,:,i,j),deg);
        end
    end

end