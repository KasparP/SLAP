function [RS, CS] = getTranslationByStrips(IM, refIM)
%align 

IM = IM-repmat(nanmean(IM,2), 1, size(IM,2), 1);
refIM = refIM - repmat(nanmean(refIM,2),1,size(refIM,2));

blocksize = 6;
res= size(IM);
nIMs = res(3);

Nstrips = ceil(res(1)/(blocksize));
xbounds = round(linspace(1, res(1), Nstrips+3));
centers = (xbounds(2:end-2) + xbounds(3:end-1))/2;
offset = nan(Nstrips, 2, nIMs); %X offsets
GoF = nan(Nstrips,nIMs);
for stripN = 1:Nstrips
    stripRef = refIM(xbounds(stripN):xbounds(stripN+3),:);
    for imN = 1:nIMs %get the rigid offset for each strip
        strip = IM(xbounds(stripN+1):xbounds(stripN+2),:, imN);
        C = normxcorr2_general(stripRef, strip, numel(strip) - (2*blocksize*size(strip,1)));
        [~,maxix] = max(C(:));
        [i,j] = ind2sub(size(C), maxix);
        i = i-ceil(size(C,1)/2); j = j-ceil(size(C,2)/2);
        
        offset(stripN,1,imN) = min(blocksize,max(-blocksize,i)); offset(stripN,2,imN) = j;
        GoF(stripN,imN) = nansum(strip(:).^2);
    end
end

RS = nan(res(1), nIMs);
CS = nan(res(1),nIMs);
for imN = 1:nIMs
    if all(all(isnan(IM(:,:,imN))))
        continue
    end
    rowshift = offset(:,1,imN);
    colshift = offset(:,2,imN);
    weight = GoF(:,imN);
    weight = weight-min(weight);
    
    meanRow = sum(rowshift.*weight)./sum(weight);
    %rowshift = medfilt2(rowshift, [3 1]);
    cs = csaps([1 centers res(1)],[-meanRow rowshift' -meanRow],1e-9,[], [nanmax(weight)/10; weight ;nanmax(weight)/10]);
    RS(:,imN) = fnval(cs, 1:res(1));
    
    meanCol = sum(colshift.*weight)./sum(weight);
    %colshift = medfilt2(colshift, [3 1]);
    cs = csaps([1 centers res(1)],[-meanCol colshift' -meanCol],1e-9,[], [nanmax(weight)/10; weight ;nanmax(weight)/10]);
    CS(:,imN) = fnval(cs, 1:res(1));
end
% 
% figure('name', 'Before Align'), imagesc(trimmean(IM, 20, 3))
% figure('name', 'After Align'), imagesc(trimmean(IM2, 20, 3))
end