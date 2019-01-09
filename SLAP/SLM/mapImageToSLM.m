function im_SLM = mapImageToSLM (im_SI, SIPixelToRefTransform, scanfield)
slmResolution = 512;
[slmGridXX,slmGridYY] = meshgrid(1:slmResolution,1:slmResolution);
[slmGridRefSpaceXX,slmGridRefSpaceYY] = xformMesh(slmGridXX,slmGridYY,SIPixelToRefTransform);

[frameGridRefSpaceXX,frameGridRefSpaceYY] = scanfield.meshgrid();

im_SLM = interp2(frameGridRefSpaceXX,frameGridRefSpaceYY,double(im_SI),slmGridRefSpaceXX,slmGridRefSpaceYY,'nearest',0);
end