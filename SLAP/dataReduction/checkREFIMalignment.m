function checkREFIMalignment (refIM, data)

%check the absolute alignment of the reference image by comparing a
%two-photon image that has been processed to a raw image of the aperture
%trransformed to galvo voltage space

%call SLAPmi_refIM on your stack imaged through the aperture

coords = refIM.M.coords;
centerZ = ceil(size(refIM.data,4)/2);

apIM = data.aperture.raw;

[Xmesh, Ymesh] = ndgrid(coords.X, coords.Y);
Xgrid = reshape(data.v2px(Xmesh(:),Ymesh(:)), size(Xmesh));
Ygrid = reshape(data.v2py(Xmesh(:),Ymesh(:)), size(Xmesh));
GI = griddedInterpolant(apIM{2});  
IM = reshape(GI(Xgrid(:), Ygrid(:)), size(Xgrid));


figure, imshow(IM,[])
figure, imshow(refIM.data(:,:,centerZ),[])
keyboard