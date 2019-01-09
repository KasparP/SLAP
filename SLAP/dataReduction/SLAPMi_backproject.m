function SLAPMi_backproject(scandata, P)
%shows a backprojection of scandata

if nargin<2
    P = linePSF_full(scandata);
end

D = [scandata.frames.pmtData];
D = nanmean(D,2);
D(isnan(D)) = 0;

BP = D'*P.P;
figure('Name', ['Backprojection for ' scandata.metadata.fileNameStem int2str(scandata.metadata.fileNameApp)]),
if length(P.coords{3})>1
     imshow3D(reshape(BP,length(P.coords{1}), length(P.coords{2}),length(P.coords{3})));
else
    imagesc(reshape(BP,length(P.coords{1}), length(P.coords{2})));
end