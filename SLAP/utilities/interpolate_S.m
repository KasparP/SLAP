function Sout = interpolate_S (S, Zcenter, nZ)
%interpolates a segmentation S onto the Z-locations defined by Zcenter and nZ (specified in pixels- requires uniform spacing)
in_sz = size(S.IM);
out_sz = [size(S.IM,1) size(S.IM,2) nZ];

planes = Zcenter-floor(nZ/2):Zcenter+floor(nZ/2);
Sout.IM = reshape(interp1(reshape(S.IM, size(S.IM,1)*size(S.IM,2), [])', planes)', size(S.IM,1), size(S.IM,2), length(planes));  %S.IM(:,:,planes)
Sout.M.coords.Z = interp1(S.M.coords.Z, planes);

if Zcenter == ceil(in_sz(3)/2);
    select = false(size(S.IM));
    select(:,:,planes) = true;
    Sout.seg = S.seg(select,:);
    return
end

[i, j, v] = find(S.seg);
[x, y, z] = ind2sub(in_sz, i);
z = z - Zcenter + (nZ+1)/2;
valid = z>=1 & z<=nZ;
x = x(valid); y = y(valid); z = z(valid); j = j(valid); v = v(valid);

%Use accumulation trick to interpolate the segmentation matrix at superresolution in Z
i_lo = sub2ind(out_sz, x, y, floor(z));
i_hi = sub2ind(out_sz, x, y, ceil(z));
i_frac = rem(z,1);
Sout.seg  =  sparse([i_hi i_lo],[j j],[i_frac.*v (1-i_frac).*v], prod(out_sz), size(S.seg,2));
end