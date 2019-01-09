function [lut, M2] = SLM_lut(mapSlm)
%makes a lookuptable from the pixelwise SLM map obtained with Georg's user
%function

M2 = nan(size(mapSlm));

%for every pixel in the SLM
lutmax = nan(size(mapSlm,1), size(mapSlm,2));
lutmin = nan(size(mapSlm,1), size(mapSlm,2));
for x = 1:size(mapSlm,1)
    x
    for y = 1:size(mapSlm,2)
        trace = smooth(squeeze(sum(sum(mapSlm(min(end,max(1,[x-1:x+1])), min(end,max(1,[y-3:y+3])), :),1),2)), 25);
        trace = trace-min(0,min(trace));
        M2(x,y,:) = trace./max(trace);
        if any(~isnan(trace))
            %find the maximum in the range 120-180
            [maxval,maxix] = max(trace(120:180));
            maxix = maxix+120-1;
            lutmax(x,y) = maxix; % + centroid(trace(maxix + [-20:20]));
            
            [minval,minix] = min(trace(180:255));
            minix = minix+180-1;
            lutmin(x,y) = minix;% + centroid(-trace(minix + [-20:20]));
        end
    end
end
lutmin = inpaint_nans(lutmin);
lutmax = inpaint_nans(lutmax);

lutmin = medfilt2(lutmin, [5 5], 'symmetric');
lutmax = medfilt2(lutmax, [5 5], 'symmetric');

lut = cat(3, lutmin, lutmax);

%update calibration file
[calibfn, calibdr] = uigetfile('E:\SLAPmidata\Calibration\calibration.cal', 'Select a filename to save calibration');
if calibfn
    calib = [];
    load([calibdr filesep calibfn], '-mat');
    calib.SLM.lut = lut;
    calib.SLM.T = ones(size(calib.SLM.lut)); %transmission in on and off state
    for x = 1:size(calib.SLM.T,1)
        for y = 1:size(calib.SLM.T,2)
            calib.SLM.T(x,y,1) = M2(x,y, calib.SLM.lut(x,y,1))./M2(x,y, calib.SLM.lut(x,y,2));
        end
    end
    calib.SLM.T(:,:,1) = medfilt2(calib.SLM.T(:,:,1), [3 3], 'symmetric');
    save([calibdr filesep calibfn], 'calib');
end
end

function peak = centroid(I)
    I = I(:) - min(I);
    ix = sum((1:length(I))'.*I)./sum(I);
    peak = ix - (length(I)+1)/2;
end