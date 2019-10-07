function [ Y, P0, P1, masks ] = phase_correction_poly( Y0, varargin )
%PHASE_CORRECTION_POLY  Background phase correction by polynomial fit
%
%   Y = PHASE_CORRECTION_POLY( Y0 );
%

%   jfpva (joshua.vanamerom@kcl.ac.uk)
%   Tom Roberts (t.roberts@kcl.ac.uk)


%% Preprocessing 

% Dimensions
[nRow,nCol,nSl,nFrame] = size( Y0 );


%% Parse Inputs

p = inputParser;

default.H           = true( size( Y0 ) );
default.U           = true( size( Y0 ) );
default.voxSize     = [1,1,1];
default.polyOrder   = 3;

addRequired(  p, 'Y0', ...
    @(x) validateattributes( x, {'numeric'}, {},mfilename) );

addParameter(  p, 'heartmask', default.H, ...
    @(x) validateattributes( x, {'logical'}, {'size',[nRow,nCol,nSl,NaN]}, mfilename));
addParameter(  p, 'uterusmask', default.U, ...
    @(x) validateattributes( x, {'logical'}, {'size',[nRow,nCol,nSl,NaN]}, mfilename));
addParameter( p, 'voxelsize',  default.voxSize, ...
    @(x) validateattributes( x, {'numeric'}, {'numel',3}, mfilename) );
addParameter( p, 'polyorder',  default.polyOrder, ...
    @(x) validateattributes( x, {'numeric'}, {'scalar','nonnegative'}, mfilename) );

parse( p, Y0, varargin{:} );

H           = p.Results.heartmask;
U           = p.Results.uterusmask;
voxSize     = p.Results.voxelsize;
polyOrder   = p.Results.polyorder;


%% Mask

% Heart Mask
if size(H,4) ~= nFrame
    H = padarray( H, [0,0,0,size(Y0,4)-size(H,4)], 'post', 'replicate' );
end

% Uterus Mask
if size(U,4) ~= nFrame
    U = padarray( U, [0,0,0,size(Y0,4)-size(U,4)], 'post', 'replicate' );
end

% Background Mask 
B = abs(Y0)<prctile(abs(Y0(U(:)&~H(:))),2);

% Combined Mask
M = U & ~H & ~B; 
M1 = M; %TAR --- store for debugging


%% Estimate Constant Background Phase

P0 = median(angle(Y0(M(:))));

% TAR --- Estimate Constant Background Phase for each slice / dynamic
% for tt = 1:size(Y0,4)
%     for ss = 1:size(Y0,3)
%         
%         P0(ss,tt) = median(angle(Y0(M(:,:,ss,tt))));
%         
%     end
% end


%% Local Temporal Variation

n = 3;
V = sqrt( 1/(n-1) * convn( abs( (Y0 - convn(Y0,ones(1,1,1,n)/n,'same') ) ).^2, ones(1,1,1,n), 'same') );  % FIXME: edges aren't calculated correctly 
V(:,:,:,1:2) = repmat(V(:,:,:,3),[1,1,1,2]);  % FIXME: for different n
V(:,:,:,end+[-1,0]) = repmat(V(:,:,:,end-2),[1,1,1,2]);
v = V(:,:,:,3:(end-2));  % FIXME: for different n
m = M(:,:,:,3:(end-2));

% Mask Non-Dynamic Voxels
N   = V < median(v(m(:)))+2.575*std(v(m(:)));  % NOTE: assuming dynamic voxels have temporal variation higher than 99% of Normal distribution

% Update Mask
M   = N & U & ~H & ~B;


% %% TAR --- View Masks
% sNum = 6;
% dNum = 50;
% imtar([U(:,:,sNum,dNum) ~H(:,:,sNum,dNum) ~B(:,:,sNum,dNum) M1(:,:,sNum,dNum) N(:,:,sNum,dNum) M(:,:,sNum,dNum)])


%% Fit Polynomial to Background Phase

[xdc,ydc,zdc] = meshgrid(voxSize(2)*(1:nCol),voxSize(1)*(1:nRow),voxSize(3)*(1:nSl));
x = repmat(xdc,[1,1,nFrame]);
y = repmat(ydc,[1,1,nFrame]);
z = repmat(zdc,[1,1,nFrame]);
model = polyfitn( [x(M(:)),y(M(:)),z(M(:))], angle(Y0(M(:))*exp(-(1i*P0))), polyOrder ); 

P1 = reshape( polyvaln( model, [xdc(:),ydc(:),zdc(:)] ), [nRow,nCol,nSl] );


%% Corrected Signal

Y = Y0 .* exp( -( 1i*(P0+P1) ) );


%% Masks

masks.combined = M;
masks.heart    = H;
masks.uterus   = U;
masks.bkgd     = B;
masks.static   = N;


end  % phase_correction_poly(...)