function [Pm, support]= applymotion(P, motion, Sz, Z)
%Applies a rigid motion transform to projection matrix P

%inputs:    
    % P:        Projection matrix, #lxls x #pixels
    % seg:      sparse image or segmentation, #pixels x #segments
    % motion:   motion vector, 1x3, in pixels
    % Z:        center plane of P on reference image in pixels, either 1 x 1  or  1 x #lxls 
    % Sz:      size of reference image in Z, in pixels
%outputs
    % Pm:       Shifted projection matrix
    % support:  if nargout==2, crop the projection matrix and return the support on S

cropping = nargout==2;

Psz = [length(P.coords{1}), length(P.coords{2}), length(P.coords{3})];
Pcenter = ceil(Psz(3)/2);
if nargin<4 || isempty(Z)
    Z = ceil(Sz/2);
end

[i, j, v] = find(P.P);
[x, y, z] = ind2sub(Psz, j);
x = round(x - motion(1));
y = round(y - motion(2));
if numel(Z==1)
    z = z - motion(3) + Z - Pcenter;
else
    z = z - motion(3) + Z(i) - Pcenter;
end

if cropping
     zmin = floor(min(z))-1;
     zmax = ceil(max(z)-zmin);
else %not cropping
    %output onto the grid of seg
     zmin = 0;
     zmax = Sz;
end

Psz(3) = Sz; %we move to the size of the sample grid
valid = x<=Psz(1) & y<=Psz(2) & z<=Psz(3) & x>=1 & y>=1 & z>=1;
x = x(valid); y = y(valid); z = z(valid); i = i(valid); v = v(valid);

j_lo = sub2ind(Psz, x, y, floor(z));
j_hi = sub2ind(Psz, x, y, ceil(z));
j_frac = rem(z,1);

%Use accumulation trick to interpolate the projection matrix at superresolution in Z
Pm.P  =  sparse([i i],[j_hi j_lo],[j_frac.*v (1-j_frac).*v], size(P.P,1), prod(Psz));
Pm.coords = P.coords; Pm.coords{3} = 1:Sz;
end