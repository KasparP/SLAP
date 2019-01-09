function R = interpolate_refIM (R, Zcenter, nZ, dX,dY)
if nargin<4 %we don't need to translate for non-masked data
    dX = 0; dY = 0;
end
%interpolates a refIM R onto the Z-locations defined by Zcenter and nZ (specified in pixels- requires uniform spacing)
planes = Zcenter + (-((nZ-1)/2):((nZ-1)/2));
%planes = Zcenter-floor(nZ/2):Zcenter+floor(nZ/2);
R.M.coords.Z = interp1(R.M.coords.Z, planes);
R.IM = reshape(interp1(reshape(double(R.IM), [], size(R.IM,3))', planes, 'linear', 0)', size(R.data,1), size(R.data,2), nZ);
R.data = reshape(interp1(reshape(double(R.data), [], size(R.data,4))', planes, 'linear', 0)', size(R.data,1), size(R.data,2), size(R.data,3), nZ);

R.labels = imtranslate(R.labels); %REIMPLEMENT THIS
keyboard
R.IM = imtranslate(R.IM, [dX dY]);
R.data = reshape(imtranslate(squeeze(R.data), [dX dY]), size(R.data));