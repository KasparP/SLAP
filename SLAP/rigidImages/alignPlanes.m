function refIM = alignPlanes(refIM)
    F2 = fft2(refIM.IM(:,:,1));
    
    rowshift = zeros(size(refIM.IM,3),1);
    colshift = rowshift;
    for z = 2:size(refIM.IM,3)
        F1 = F2;
        F2 = fft2(refIM.IM(:,:,z));
        output = dftregistration(F1,F2,4);
        
        rowshift(z) = output(3);
        colshift(z) = output(4);
    end
    
    rowshift = cumsum(rowshift);
    colshift = cumsum(colshift);
    rowshift = rowshift - round(4*mean(rowshift))/4;
    colshift = colshift - round(4*mean(colshift))/4;
    
    %IM2 = zeros(size(refIM.IM));
    for z = 1:size(refIM.IM,3)
       refIM.IM(:,:,z) = imtranslate(refIM.IM(:,:,z), [colshift(z) rowshift(z)]);
       for ch = 1:size(refIM.data,4)
           refIM.data(:,:,z,ch)= imtranslate(refIM.data(:,:,z,ch), [colshift(z) rowshift(z)]);
       end
    end
end