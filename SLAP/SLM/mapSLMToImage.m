function im_SI = mapSLMToImage (im_SLM, SIPixelToRefTransform, scanfield)
slmResolution = 512;
[slmGridXX,slmGridYY] = meshgrid(1:slmResolution,1:slmResolution);
[slmGridRefSpaceXX,slmGridRefSpaceYY] = xformMesh(slmGridXX,slmGridYY,SIPixelToRefTransform);

[frameGridRefSpaceXX,frameGridRefSpaceYY] = scanfield.meshgrid();

Iplnt = scatteredInterpolant(slmGridRefSpaceXX(:),slmGridRefSpaceYY(:),double(im_SLM(:)),'nearest');
im_SI = Iplnt(frameGridRefSpaceXX,frameGridRefSpaceYY);
%im_SI = interp2(slmGridRefSpaceXX,slmGridRefSpaceYY,double(im_SLM),frameGridRefSpaceXX,frameGridRefSpaceYY,'nearest',NaN);
end
