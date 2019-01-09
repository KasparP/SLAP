function SLAPMi_alignTwoREFChannels
basedir = 'E:\SLAPMidata\';

[fnR, drR] = uigetfile([basedir '*.mat'], 'select your Red channel REFim');
[fnG, drG] = uigetfile([drR '*.mat'], 'select your Green channel REFim');


SR = load([drR fnR]);
SG = load([drG fnG]);

R = squeeze(SR.refIM.data);
G = squeeze(SG.refIM.data);


RXY = max(R,[],3);
GXY = max(G,[],3);


%align max intensity projection in 2D
%[rshift, cshift] = registerXcorr(R2,G2, 30);
[Xshift1, Yshift1] = registerDFT(RXY,GXY);

%shift the red image
R2 = imtranslate(R, -[Yshift1 Xshift1 0]);

RXZ = squeeze(max(R2,[],2));
RYZ = squeeze(max(R2,[],1));

GXZ = squeeze(max(G,[],2));
GYZ = squeeze(max(G,[],1));

[Xshift2, Zshift1] = registerDFT(RXZ,GXZ);
[Yshift2, Zshift2] = registerDFT(RYZ,GYZ);

Rfinal = imtranslate(R, [-(Yshift1+Yshift2) -(Xshift1+Xshift2) (Zshift1+Zshift2)/2]);

%save the merged REFIM
refIM = SG.refIM;
refIM.data = repmat(refIM.data, [1 1 3 1]);
refIM.data(:,:,1,:) = Rfinal;
refIM.data = im2uint8(refIM.data);

save([drG fnG(1:end-4) 'MERGED'], 'refIM', '-v7.3'); %mat file)
options.overwrite = true; options.color = true;
errorcode = saveastiff(refIM.data, [drG fnG(1:end-4) 'MERGED.tif'], options); %tiff file for easy reading
if errorcode
    keyboard
end
disp('Done merging!')
end


    function [rshift, cshift] = registerXcorr(A,B, maxshift)
        requiredNumberOfOverlapPixels = prod(max(1,min(size(A), size(B))-maxshift));
        cc = normxcorr2_general(A,B,requiredNumberOfOverlapPixels);
        [rshift,cshift] = find(cc==max(cc(:)));
        rshift = rshift-size(A,1);
        cshift = cshift-size(A,2);
    end
    
    function [rshift, cshift] = registerDFT(A,B)
        output = dftregistration(fft2(A), fft2(B),2);
        rshift = output(3);
        cshift = output(4);
    end
    