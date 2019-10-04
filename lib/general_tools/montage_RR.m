function [ bigImage ] = montage_RR( input_args,  Username, Userlims)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

A = size(input_args);
if length(A) == 3
    B = reshape(input_args, A(1), A(2), 1, A(3));
    C = input_args(:,:,1);
elseif length(A) == 4
    B = reshape(input_args, A(1), A(2), 1, A(3)*A(4));
    C = input_args(:,:,1,1);
% TAR -- accept 2D data, ie, one slice
elseif length(A) == 2
    B = input_args;
    C = input_args;
end

if nargin == 1
    Username = '';
    Usermax = max(input_args(find(input_args)));
%     Usermax = max(C(find(C)));
    Userlims = [0 0.8*Usermax];
elseif nargin ==2
    if isnumeric(Username)
        Userlims = Username;
        Username = '';
    else
        Usermax = max(C(find(C)));
        Userlims = [0 0.8*Usermax];
    end
elseif nargin == 3
%     Usermax = max(C(find(C)));
%     Userlims = [0 Usermax];
end


figure('Name', Username);

h = montage(B);
set(gcf, 'OuterPosition', [50 50 1000 1000])



set(gca, 'CLim', Userlims)
colormap('gray');

hh = get(h);
bigImage = hh.CData;

end

