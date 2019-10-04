function [Vworld, Vxyz] = gradfirstmom_mps2world( Vmps, gc )

%% GRADFIRSTMOM_MPS2WORLD	convert gradient moment vectors from scanner to world coordinates
%
%   GRADFIRSTMOM_MPS2WORLD( Vmps, gc )
%       Convert gradient first moment vector from scanner (MPS) to world
%       coordinates (and XYZ coordinates)
%       - uses the O_MATRIX from the Philips goalc code associated with the 
%         stack to transform from scanner to world coordinates
%
% Input:
%       - Vmps          First moments measured in GVE
%       - gc            goalc structure containing m_/p_/s_/orient
%
% Output:
%       - Vworld        First moments in world coordinates
%       - Vxyz          First moments in scanner (xyz) coordinates
%
% TODO:
% - eventually change so that it can calculate Vmps from goalc structure
%
%

% Tom Roberts (t.roberts@kcl.ac.uk)

%% First moments
Vm = Vmps(1);
Vp = Vmps(2);
Vs = Vmps(3);

Vmps = [Vm; Vp; Vs];


%% F --- First moment orientations
% think I can use this to automatically determine gradient moment directions
F.m = [gc.mc0_str, gc.m0_str ];  % should be [-ve, +ve]
F.s = [gc.s_ex_str, gc.r0_str ]; % should be [+ve, -ve]

if sign(F.m) ~= [-1, 1]
    warning('Readout gradients are orientated in an expected way!');
end

if sign(F.s) ~= [1, -1]
    warning('Slice-select gradients are orientated in an unexpected way!');

    % update Vmps if slice-select is flipped
    if sign(F.s) == [-1, 1]
        Vs = Vs * -1;
        Vmps = [Vm; Vp; Vs];
        warning('Slice-select gradient flipped: Vs = -Vs');
    end
end


%% O_MATRIX --- converts direction of slice from MPS to xyz coordinates
%- get from goalc: 'O_MATRIX `locations[0]'
%- eg: for transverse slice:
%- O_MATRIX = [m_orient; p_orient; s_orient];
%- m_orient = [0 -1 0]; <- means M in -y direction (= LR in Philips)
%- p_orient = [1 0 0];  <- means P in +x direction (= PA in Philips)
%- s_orient = [0 0 -1]; <- means S in -z direction (= HF in Philips)
m_orient = gc.m_orient;
p_orient = gc.p_orient;
s_orient = gc.s_orient;
O = [m_orient; p_orient; s_orient];


%% Define Gradient Moment direction
%- relationship between direction of gradients and +ve/-ve velocity encoding
%- Either positive or negative
%- For affine approach: worked out that needs to be negative by matching to equivalent QFlow
% D = [1, 0, 0; 0 1, 0; 0 0 1];
D = [-1, 0, 0; 0 -1, 0; 0 0 -1];


%% First moments in xyz coordinates
Vxyz = D * O' * Vmps;


%% Transform from scanner to world
%- Philips defines scanner as:
%- +x = PA
%- +y = RL
%- +z = FH (IS)

%- Sensible World/Patient definition (ie: RAI):
%%%%% I think this might need to match .nii coord system? ie: LAS = radiological?
%%%%% nb: when I corrected Josh's .nii output, it puts .nii into neurological order = RAS
%- +x` = RL
%- +y` = AP
%- +z` = FH (IS)

%- conversion from Philips to normal world:
%- scanner : world
%-       x : -y`
%-       y :  x`
%-       z :  z`

% think these Cprimes are not incorrectly defined... 
% despite "RAI" below working in affine data
% First Cprime transform actually gives LPI (rather than RAI)
Cprime = [0 -1 0; 1 0 0; 0 0 1];  %<- RAI --- this works with affine approach
% Cprime = [0 -1 0; 1 0 0; 0 0 -1]; %<- RAS
% Cprime = [0 1 0; 1 0 0; 0 0 -1];  %<- LAS --- matches Josh's default .nii output

% correct transforms
% Cprime = [0 1 0; -1 0 0; 0 0 -1];  %<- LAS --- matches Josh's default .nii output

%% Perform transform to world coordinates
%- currently, based on Cprime above:
%- Vworld(1) = RL
%- Vworld(2) = AP
%- Vworld(3) = FH
Vworld = Cprime * Vxyz;


% end gradfirstmom_mps2world(...).m
end