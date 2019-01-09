function figure_uncaginglxlDFFplot(scandata)

tau = 100;
stimLength = 10;
[stimArtifact, stimStart] = estimateStimArtifact(scandata,stimLength);

y = [scandata.frames.pmtData];
yFix = y-stimArtifact;

[F, dFF] = fastDFF (yFix, tau);

dFF_0 = max(0, dFF - dFF(:,stimStart));
figure, imagesc(dFF_0(:,stimStart-200+1:stimStart+1000),[prctile(dFF_0(:),10),prctile(dFF_0(:),99.9)]);
%%

c = PaperColorMap(); 
% colormap(jet); colorbar; 
colormap(c); colorbar

%%
%draw rectangles
hold on, 
edges = [0 250 500 750 1000];
C = [0 1 0; 0.2 0.2 1; 1 0 0; 0.8 0.8 0.2];
for l = 1:4
    line([0 size(dFF,2) size(dFF,2) 0 0], [edges(l) edges(l) edges(l+1) edges(l+1) edges(l)], 'color', C(l,:))
end
%draw stimulus time
line([200 200], [0 size(dFF_0,1)], 'linestyle', ':', 'color', 'r', 'linewidth', 2)

end