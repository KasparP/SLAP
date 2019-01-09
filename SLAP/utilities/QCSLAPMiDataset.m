function [dataset, selection] = QCSLAPMiDataset (dataset)
%select data that passes QC
selected = 100*ones(1,size(dataset.dPhotons,3));
notselected = zeros(1,size(dataset.dPhotons,3));
BDpt = [];
colors = hsv(8);
hF = figure('Name', 'Select trials to include');
hAx = axes('parent', hF);
hIm =imagesc(squeeze(nansum(dataset.dPhotons,1))');
hold on;
sInclude = scatter(100*ones(1,length(selected)), 1:length(selected), 'markerFaceColor', 'g', 'markeredgecolor', 'g');
sExclude = scatter(nan(1,length(selected)), 1:length(selected), 'markerFaceColor', 'r', 'markeredgecolor', 'r');
for stimIX = 1:8
    S = dataset.stimulus.stim==stimIX;
    sStim(stimIX) = scatter((size(dataset.dPhotons,2)-100)*ones(1,sum(S)), find(S), 'markerFaceColor', colors(stimIX,:), 'markeredgecolor', colors(stimIX,:));
end
xlabel('Time')
ylabel('Trials')

set(hF, 'WindowButtonDownFcn', @selectBDF)
% set(hF, 'WindowButtonMotionFcn', @selectBMF)
% set(hF, 'WindowButtonUpFcn', @selectBUF)

waitfor(hF);
selection = selected>0;

if nargout==1
    dataset.stimTime = dataset.stimTime(selection);
    dataset.dPhotons = dataset.dPhotons(:,:,selection);
    dataset.F0_photons = dataset.F0_photons(:,selection);
    dataset.X = dataset.X(:,:,selection);
    dataset.Z = dataset.Z(selection);
    dataset.Zcenter = mode(dataset.Z) + ceil(size(dataset.refIM.IM,3)/2);
    dataset.spikes = dataset.spikes(:,:,selection);
    dataset.stimulus.stim = dataset.stimulus.stim(selection);
    dataset.stimulus.stimTime = dataset.stimulus.stimTime(selection);
    dataset.stimulus.timeError = dataset.stimulus.timeError(selection);
    dataset.filenames = dataset.filenames(selection);
else
    dataset = [];
end
    function selectBDF(evnt,data)
        cp = get(hAx,'CurrentPoint');
        cp = cp(1, 1:2);
        BDpt = round(cp(2));
        if BDpt>=1 && BDpt<=length(selected)
        selected(BDpt) = 100- selected(BDpt);
        notselected(BDpt) = 100- notselected(BDpt);
        set(sInclude, 'xdata', selected);
        set(sExclude, 'xdata', notselected);
        end
    end
end