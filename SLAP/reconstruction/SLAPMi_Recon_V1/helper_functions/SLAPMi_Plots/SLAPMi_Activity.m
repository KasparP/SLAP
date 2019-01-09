function SLAPMi_Activity(sys_recon)

ref2D = sum(sys_recon.input.ref_image,3);
ysz = size(ref2D);
S = 0;
for plane = 1:sys_recon.opts.Psize
S = S+ sys_recon.input.S((1:numel(ref2D)) + (plane-1)*numel(ref2D),:);
end

% colors = hsv(nStim);
hF = figure('color', 'w');
hIm = imagesc(ref2D);
hAxIm = get(hIm, 'parent'); %set(hAxIm, 'pos', [0.1 01 1 1]);
% hold(hAxIm,'on');

hF2 = figure;
hAxActivity= axes('Pos', [0.1 0.1 0.8 0.8], 'tickdir', 'out', 'linewidth', 1.5);
set (hF, 'WindowButtonDownFcn', @bdf);

    function bdf(src, event)
        cp = get(hAxIm,'CurrentPoint');
        cp = cp(1, 1:2);
        %If the point is within the axes
        if all(cp>=1 & cp<=ysz)
            ptPos = round([cp(1) cp(2)]);
            ind = sub2ind(size(ref2D),ptPos(2), ptPos(1));
            bestSeg = find(S(ind,:)>0);
            if ~isempty(bestSeg)

                set(hF,'name', ['Segment:' int2str(bestSeg)]),
                set(hF2,'name', ['Segment:' int2str(bestSeg)]),
                
                %Activity plot
                cla(hAxActivity);
                plot(hAxActivity, sys_recon.output.F(bestSeg,:)');
            end
        end
    end


end
