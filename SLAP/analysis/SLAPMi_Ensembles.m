function SLAPMi_Ensembles(dataset, opts)

dr = 'W:\ForPublication\InVivoGluSnFR\SLAPMiTuningDatasets';
[opts.fn, opts.dr] = uiputfile([dr filesep dataset.filenames{1}(1:20) '.mat']);

tau = dataset.solverOpts.tau;
theta = [1 -exp(-1/tau)];

prePeriod = 300:1000;
stimPeriod = 1100:3000;
postPeriod = 3300:3949;

if nargin<2 || ~isfield(opts, 'selection')
    [~,selection] = QCSLAPMiDataset(dataset);
elseif all(opts.selection==1)
    selection = true(size(dataset.stimulus.stim));
else
    selection = opts.selection;
end

%save memory
if isfield(dataset.refIM, 'data')
    dataset.refIM = rmfield(dataset.refIM, 'data');
    dataset = rmfield(dataset, {'spikes', 'P'});
end

%clipping to first XXX repetitions
triSel = find(selection);
dataset.stimulus.stim = dataset.stimulus.stim(triSel);
dataset.stimulus.timeError = dataset.stimulus.timeError(triSel);
dataset.stimulus.stimTime = dataset.stimulus.stimTime(triSel);
dataset.filenames = dataset.filenames(triSel);
dataset.dPhotons = dataset.dPhotons(:,:,triSel);

dataset.dPhotonSpikes = filter(theta,1,dataset.dPhotons,[],2);
dataset.Ysum = dataset.Ysum(:,triSel);
dataset.X = dataset.X(:,:,triSel);
dataset.Z = dataset.Z(triSel);
dataset.stimTime = dataset.stimTime(triSel);
% dataset.spikes = dataset.spikes(:,:, triSel);
dataset.alignError = dataset.alignError(:,triSel);
dataset.motion.shifts = dataset.motion.shifts(:,:,triSel);
dataset.motion.error = dataset.motion.error(:,triSel);

%select ROIs for analysis
%minimum response
meanAmp = nanmean(nanmean(dataset.dPhotons(:,stimPeriod, :),2),3) - nanmean(nanmean(dataset.dPhotons(:,prePeriod, :),2),3);
ampThresh = 0.2;
stims = unique(dataset.stimulus.stim);
nStim = length(stims);
for stimix = nStim:-1:1
    Xn(:,stimix) = sum(isnan(dataset.X(:,1,dataset.stimulus.stim==stims(stimix))),3);
end
Xn = max(Xn,[],2);
segSel = meanAmp>ampThresh  & Xn<3;
nMin = 200;
if sum(segSel)<nMin 
    %add in segments by total activity
    meanAmp2 = nanmean(nanmean(dataset.dPhotons(:,stimPeriod, :),2),3);
    ampThresh2 = prctile(meanAmp2(~isnan(meanAmp2)), 100*min(1, (sum(~isnan(meanAmp2))-nMin)/sum(~isnan(meanAmp2))));
    segSel = segSel | (meanAmp2>ampThresh2 & Xn<3); 
end
segSel = segSel(1:end-2);
nSeedSel = sum(segSel);
%visualize selected segments
colorSegs(dataset,segSel, ones(sum(segSel),1));

dataset.Ysum = squeeze(nansum(dataset.dPhotons,1));

%apply selection
dataset.dPhotons = dataset.dPhotons(segSel,:,:);
dataset.X = dataset.X(segSel,:,:);
%dataset.spikes = dataset.spikes(segSel,:, :);

%bleaching correction
resp = dataset.dPhotons;
normfactor =sqrt(nanmean(resp(:,50:end-50,:).^2, 2)); normfactor = normfactor./nanmedian(normfactor,3); %correct for bleaching across trials while maintaining rough photon counts
resp = resp./normfactor;
[summary.d,summary.dX,summary.dY] = find_distances(dataset, segSel);

%%%%%%%%%%%%  ANALYSES:

summary.CCevoked = figure_crossCorrelations(resp(:,stimPeriod,:), summary.d, 'evoked', summary.dX, summary.dY, dataset.stimulus.stim);
summary.CCspont = figure_crossCorrelations(resp(:,prePeriod,:), summary.d, 'spontaneous',  summary.dX, summary.dY, dataset.stimulus.stim);

summary.NMFglobal = figure_NMF_global(dataset, resp, stimPeriod, prePeriod, postPeriod, segSel,summary.d);

summary.Gshuff =  globalModesTrialShuffled(dataset, resp, stimPeriod, prePeriod, postPeriod, [], segSel);

summary.ensembles = figure_ensembles(dataset, HPfilter(resp), stimPeriod, prePeriod, postPeriod);


summary.PCAperStim = figure_PCA_perStim(dataset, resp, stimPeriod, prePeriod, postPeriod, true); %angles between first PCA components for each stimulus direction, and spontaneous 


summary.freq = figure_frequencyAnalysis(dataset, summary.ensembles.e,summary.ensembles.eShuff, summary.Gshuff,prePeriod, stimPeriod, meanAmp(segSel));

summary.PCAglobal = figure_PCA_global(dataset, resp, stimPeriod, prePeriod, postPeriod, true, segSel);  %

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

summary.filenames = dataset.filenames;
summary.dr = dataset.dr;
summary.stim = dataset.stimulus.stim;

save([opts.dr filesep opts.fn], 'summary', '-v7.3');
end

function Gshuff = globalModesTrialShuffled(dataset, resp, stimPeriod, prePeriod, postPeriod, ORI, segSel)
    repeats = 10;
    nSVD = 10;
    D = resp(:,50:end-50,:);
    discardTrials =squeeze(max(sum(isnan(D),1),[],2))>3;
    Dtmp = D(:,:,~discardTrials);
    discardSegs = any(any(isnan(Dtmp),2),3); clear Dtmp;
    D = D(~discardSegs,:,:);
    Gshuff.selectedSegs = ~discardSegs;
    discardTrials =squeeze(max(sum(isnan(D),1),[],2))>0;
    D = D(:,:,~discardTrials);
    Gshuff.selectedTrials = ~discardTrials;
    Gshuff.stim = dataset.stimulus.stim(~discardTrials);
    Gshuff.D = D;

    %shuffle stimuli to get a null for the amount of variance explained by
    %stimulus identity
    mStim = nan(size(D,1),size(D,2), 8);
    for repeat = 1:repeats
        Dsub = D;
        for m = 1:size(Dsub,3)
            Dsub(:,:,m) = Dsub(:, randperm(size(Dsub,2)), m);
        end
        Dsub = Dsub(:,:,randperm(size(Dsub,3)));
        for stim = 1:8
            sel = find(Gshuff.stim==stim);
            mStim(:,:,stim) = nanmean(Dsub(:,:,sel),3) - nanmean(nanmean(Dsub(:,[prePeriod]-49,sel),3),2);
            Dsub(:,:,sel) = Dsub(:,:,sel) - mStim(:,:,stim);
        end
        Gshuff.varDsub(:,repeat) = var(reshape(Dsub, size(Dsub,1),[]),[],2);
    end
    
    %Stimulus subtraction
    Dsub = D;
    mStim = nan(size(D,1),size(D,2), 8);
    for stim = 1:8
        sel = find(Gshuff.stim==stim);
        mStim(:,:,stim) = nanmean(D(:,:,sel),3) - nanmean(nanmean(D(:,[prePeriod]-49,sel),3),2);
        Dsub(:,:,sel) = Dsub(:,:,sel) - mStim(:,:,stim);
    end
    D = D-mean(D,2);
    Dsub = Dsub-mean(Dsub,2);
    
    %compute % of variance explained by stimulus identity, versus % of
    %variance explained by global mode
    Gshuff.varTotal = var(reshape(D,size(D,1),[]), [],2);
    Gshuff.varDsub_obs = var(reshape(Dsub,size(D,1),[]), [],2);
   
    Gshuff.G = nan(repeats,size(Dsub,2),size(Dsub,3));
    Gshuff.varGsub = nan(size(D,1), repeats);
    for repeat = 0:repeats
        disp(['Repeat #: ' int2str(repeat)]);
        %randomize trials (within stimulus??)
        if repeat>0
            for seg = 1:size(Dsub,1)
                for stim = 1:8
                    stimSel = find(Gshuff.stim==stim);
                    p = randperm(length(stimSel));
                    Dsub(seg,:,stimSel) =  Dsub(seg,:,stimSel(p));
                    D(seg,:,stimSel) =  D(seg,:,stimSel(p));
                end
            end
        end
        
        %calculate global mode without stimulus subtraction
        [U,S,V] = svd(reshape(D,size(D,1),[]),'econ');
        for comp = nSVD:-1:1
            corrSqToSumD(comp) = corr(V(:,comp), sum(reshape(D,size(D,1),[]),1)').^2;
        end
        [~,gComp] = max(corrSqToSumD);
        Gsub = reshape(D,size(D,1),[]) - (U(:,gComp)*S(gComp,gComp)*V(:,gComp)');
        if repeat==0
            Gshuff.varGsub_obs = var(Gsub, [],2);
            Gshuff.pred = reshape(U(:,gComp)*S(gComp,gComp)*V(:,gComp)', size(Gshuff.D));
        else
            Gshuff.varGsub(:,repeat) = var(Gsub, [],2);
        end
        
        %calculate global mode with stimulus subtraction
        [U,S,V] = svd(reshape(Dsub,size(Dsub,1),[]),'econ');
        w= U(:,1:nSVD)./sign(mean(U(:,1:nSVD),1));
        s= diag(S(1:nSVD,1:nSVD))./sqrt(size(V,1));
        h = V(:,1:nSVD)./sign(mean(U(:,1:nSVD),1));
        for comp = 1:nSVD
            corrSqToSumD(comp) = corr(V(:,comp), sum(reshape(D,size(D,1),[]),1)').^2;
        end
        [~,gComp] = max(corrSqToSumD);
        %generate the global predictor following stimulus subtraction
        if repeat==0
            Gshuff.gComp = gComp; Gshuff.S_obs = diag(s); 
            disp(['The global mode is SVD component ' int2str(gComp)]);
            for comp = 1:nSVD
                Gshuff.G_obs(:,:,comp) = reshape(h(:,comp).*s(comp), size(Dsub,2),size(Dsub,3));
                Gshuff.W_obs(:,:,comp) = w(:,comp).*s(comp);
            end
        else
            Gshuff.G(repeat,:,:) = reshape(h(:,gComp).*s(gComp), size(Gshuff.G(repeat,:,:)));
            Gshuff.W(:,repeat) = w(:,gComp).*s(gComp);
        end
    end
    clear U S V
    
    %distribution of global mode vs stimulus variance explained, per segment
    fracS_all = 1- (Gshuff.varDsub_obs./Gshuff.varTotal);
    n = numel(D); 
    p = size(D,1).*size(D,2).*(size(mStim,3)); fracS_all = fracS_all - (1-fracS_all).*(p/(n-p-1));%adjust R2
    fracG_all = 1- (Gshuff.varGsub_obs./Gshuff.varTotal);
    n = numel(D); p = (size(D,1)+size(D,2)).*(size(D,3)); fracG_all = fracG_all - (1-fracG_all).*(p/(n-p-1));%adjust R2
    
    for repeat = repeats:-1:1
        fracGshuff_all(:,repeat) = (Gshuff.varTotal-Gshuff.varGsub(:,repeat))./Gshuff.varTotal;
        n = numel(D); 
        p = (size(D,1)+size(D,2)).*(size(D,3)); fracGshuff_all(:,repeat) = fracGshuff_all(:,repeat) - (1-fracGshuff_all(:,repeat)).*(p/(n-p-1));%adjust R2
        fracSshuff_all(:,repeat) = (Gshuff.varTotal-Gshuff.varDsub(:,repeat))./Gshuff.varTotal;
        n = numel(D); 
        p = size(D,1).*size(D,2).*(size(mStim,3)); fracSshuff_all(:,repeat) = fracSshuff_all(:,repeat) - (1-fracSshuff_all(:,repeat)).*(p/(n-p-1));%adjust R2
    end
    %subtract shuffle
    Gshuff.fracG_all = fracG_all-mean(fracGshuff_all,2);
    Gshuff.fracS_all = fracS_all-mean(fracSshuff_all(:));
    Gshuff.fracSshuff_all = fracSshuff_all-mean(fracSshuff_all(:));
    figure, scatter(Gshuff.fracS_all, Gshuff.fracG_all); xlabel('Fraction of variance explained by mean stimulus response (Adjusted R^2)'); ylabel('Fraction of variance explained by top SVD component (Adjusted R^2)');
    figure, hist(Gshuff.fracG_all, 50); xlabel('Fraction of variance explained by top SVD component (Adjusted R^2)'); ylabel('Counts')
    figure, hist(Gshuff.fracS_all, 50); xlabel('Fraction of variance explained by mean stimulus response (Adjusted R^2)'); ylabel('Counts')
    
    figure, scatter(Gshuff.fracS_all, Gshuff.fracG_all)
    hold on, scatter(Gshuff.fracS_all, 0*Gshuff.fracG_all); %, 'markeredgecolor', [0.5 0.5 0.5])
    hold on, scatter(mean(Gshuff.fracSshuff_all,2), Gshuff.fracG_all); %, 'markeredgecolor', [0.5 0.5 0.5])
    xlabel('Fraction of variance explained by stimulus (Adjusted R^2)');
    ylabel('Fraction of variance explained by top SVD component (Adjusted R^2)');
    legend('Observed', 'Trials shuffled for each stimulus independently across segments', 'Timepoints shuffled across stimuli');
    
    %fraction of total variance explained (compute p value? p<<1e-5)
    Gshuff.fracS = sum(Gshuff.fracS_all.*Gshuff.varTotal)./sum(Gshuff.varTotal);
    Gshuff.fracG = sum(Gshuff.fracG_all.*Gshuff.varTotal)./sum(Gshuff.varTotal);
    
%     %TODO: tuning depth before and after global mode subtraction
%     TCD = nan(size(Gshuff.D,1), 8); TCP = TCD;
%     for stim = 1:8
%         sel = find(Gshuff.stim==stim);
%         TCD(:,stim) = mean(mean(Gshuff.D(:,stimPeriod-49,sel), 2),3) - mean(mean(Gshuff.D(:,prePeriod-49,sel), 2),3);
%         TCP(:,stim) = mean(mean(Gshuff.pred(:,stimPeriod-49,sel), 2),3) - mean(mean(Gshuff.pred(:,prePeriod-49,sel), 2),3);
%     end
%     TCsub = TCD - TCP;
%     rad = linspace(0,2*pi,8/2+1);
%     rad = repmat(rad([1:4 1:4]),size(Gshuff.D,1),1);
%     [H_d,S_d,V_d] = circular_mean(rad,TCD);
%     [H_p,S_p,V_p] = circular_mean(rad,TCP);
%     [H_sub,S_sub,V_sub] = circular_mean(rad,TCsub);
%     figure('name', 'ORI index before global mode subtraction')
%     hist(S_d,0.025:0.05:1), xlabel('ORI index'), ylabel('counts')
%     set(gca, 'ylim', [0 120])
%     h = findobj(gca,'Type','patch'); h.FaceColor = [0 0 0.5];
%     figure('name', 'ORI index after global mode subtraction')
%     hist(S_sub, 0.025:0.05:1), xlabel('ORI index'), ylabel('counts')
%     set(gca, 'ylim', [0 120])
%     h = findobj(gca,'Type','patch'); h.FaceColor = [0 0.5 1];
%     center = H_p(1);
%     Ha = H_d-center; Ha(Ha>0.5) = Ha(Ha>0.5)-1;  Ha(Ha<-0.5) = Ha(Ha<-0.5)+1;
%     Hb = H_sub-center; Hb(Hb>0.5) = Hb(Hb>0.5)-1; Hb(Hb<-0.5) = Hb(Hb<-0.5)+1;
%     figure, scatter(Ha,Hb); xlabel('preferred orientation before Global mode subtraction'); ylabel('preferred orientation after global mode subtraction');
    
    %color an image according to variance explained by global(red) and stimulus(green)
%     colors = zeros(sum(segSel),2); colors(Gshuff.selectedSegs,1) = Gshuff.fracG_all; colors(Gshuff.selectedSegs,2) = Gshuff.fracS_all; 
%     colorSegs(dataset, segSel, colors);set(gcf, 'name', 'Red: global modulation; Green: Stimulus modulation')

    %tuning of global activity predictors
    for stim = 1:8
        sel = Gshuff.stim==stim;
        TCdata = squeeze(var(Gshuff.G_obs(stimPeriod-49,sel,:),[],1) - var(Gshuff.G_obs(prePeriod-49,sel,:),[],1));  %- mean(Gshuff.G_obs(prePeriod-49,sel),1);
        Gshuff.TC_obs(stim,:) = mean(TCdata,1); Gshuff.TC_obsE(stim,:) = std(TCdata,0,1)./sqrt(size(TCdata,1));
        TCdata = squeeze(var(Gshuff.G(:,stimPeriod-49,sel),[], 2)) - squeeze(var(Gshuff.G(:,prePeriod-49,sel),[], 2)); %squeeze(mean(Gshuff.G(:,stimPeriod-49,sel),2) - mean(Gshuff.G(:,prePeriod-49,sel),2));
        Gshuff.TC(:,stim) = mean(TCdata,2);  Gshuff.TCE(:,stim) = std(TCdata,0,2)./sqrt(size(TCdata,2)); 
    end
    stim = 9;
    TCdata = squeeze(var(Gshuff.G_obs([prePeriod postPeriod]-49,:,:),[],1));
    Gshuff.TC_obs(stim,:)= mean(TCdata,1); Gshuff.TC_obsE(stim,:) = std(TCdata,0,1)./sqrt(size(TCdata,1));
    TCdata = squeeze(var(Gshuff.G(:,[prePeriod postPeriod]-49,sel),[], 2)); %squeeze(mean(Gshuff.G(:,stimPeriod-49,sel),2) - mean(Gshuff.G(:,prePeriod-49,sel),2));
    Gshuff.TC(:,stim) = mean(TCdata,2);  Gshuff.TCE(:,stim) = std(TCdata,0,2)./sqrt(size(TCdata,2)); 
    
    %FIGURE; here the null is the trial-shuffled-per-stimulus power in the global mode
    figure('name', 'Fraction increase in the power of the global mode during the stimulus period, minus Null (controls for change in brightness)'), 
    modecolors = parula(nSVD);
    for mode = 1:nSVD % Gshuff.gComp
        errorbar((Gshuff.TC_obs(1:8,mode)'-mean(Gshuff.TC(:,1:8),1))./Gshuff.TC_obs(9,mode), Gshuff.TC_obsE(1:8,mode)./Gshuff.TC_obs(9,mode), 'color', modecolors(mode,:), 'linewidth', 1.5)
        hold on
        colors = hsv(8);
        for stimix=1:8
            scatter(stimix, (Gshuff.TC_obs(stimix,mode)-mean(Gshuff.TC(:,stimix),1))./Gshuff.TC_obs(9,mode), 80, 'markerfacecolor', colors(stimix,:), 'markeredgecolor', modecolors(mode,:), 'linewidth', 1.5)
        end
        hold on, plot([0 9], [0 0], 'k:');
        xlabel('stimulus'); ylabel('\Delta power (normalized)')
        set(gca, 'xlim', [0 9], 'tickdir', 'out')
    end
    
end

function PSD = figure_frequencyAnalysis(dataset, e,eShuff, Gshuff, prePeriod, stimPeriod, amp)
tau = dataset.solverOpts.tau;
theta = [1 -exp(-1/tau)];

%frequency spectra of the top 5 brightest spines
[sorted, sortorder] = sort(amp, 'descend');
stimuli =  dataset.stimulus.stim;
% NS =  10;
% for ix = 1:NS
%     [PSD.brightSegs.total(:,ix), PSD.brightSegs.freqs] = makePSD(squeeze(dataset.dPhotons(sortorder(ix),50:end-50,:)), stimuli);
%     PSD.brightSegs.evoked(:,ix) = makePSD(squeeze(dataset.dPhotons(sortorder(ix),stimPeriod,:)), stimuli);
%     PSD.brightSegs.spont(:,ix) = makePSD(squeeze(dataset.dPhotons(sortorder(ix),prePeriod,:)), stimuli);
% end
% figure('name', ['Mean power spectra of ' int2str(NS) ' brightest segments'])
% plot(PSD.brightSegs.freqs,mean(PSD.brightSegs.total,2), 'color', 'b'); hold on, 
% plot(PSD.brightSegs.freqs,mean(PSD.brightSegs.evoked,2)); hold on, plot(PSD.brightSegs.freqs,mean(PSD.brightSegs.spont,2));
% legend({'total','evoked', 'spontaneous'})
% set(gca, 'yscale', 'log', 'xscale', 'log'); 
% figure('name', 'evoked power/spontaneous power in brightest segments'), 
% plot(PSD.brightSegs.freqs,PSD.brightSegs.evoked./PSD.brightSegs.spont); xlabel('frequency'), ylabel('power ratio')

%frequency spectrum of global events
dtype = {'ensemble', 'pca'};
for d = 1:2
    if d==1
        data = e; datanull = eShuff; stimuli = dataset.stimulus.stim;
    else
        data = Gshuff.G_obs(:,:,Gshuff.gComp); datanull = permute(Gshuff.G, [2 3 1]); stimuli = dataset.stimulus.stim(Gshuff.selectedTrials);
    end
    %deconvolve
    data = filter(theta,1,data,[],1);
    datanull =  filter(theta,1,datanull,[],1);
    
    [PSD.(dtype{d}).total, PSD.(dtype{d}).freqs{1}] = makePSD(data(1:size(e,1),:), stimuli);
    [PSD.(dtype{d}).evoked, PSD.(dtype{d}).freqs{2}] = makePSD(data(stimPeriod-49,:), stimuli);
    [PSD.(dtype{d}).spont, PSD.(dtype{d}).freqs{3}] = makePSD(data(prePeriod-49,:), stimuli);
    for repeat = 1:size(datanull,3)
        PSD.(dtype{d}).null.total(:,repeat) = makePSD(datanull(1:size(e,1),:,repeat), stimuli);
        PSD.(dtype{d}).null.evoked(:,repeat) = makePSD(datanull(stimPeriod-49,:,repeat), stimuli);
        PSD.(dtype{d}).null.spont(:,repeat) = makePSD(datanull(prePeriod-49,:,repeat), stimuli);
    end
    
    selections = {'total', 'evoked', 'spont'};
    colors = {'r','g','b'};
    for s = 1:3
        figure('name', ['Power spectra:' dtype{d} ', ' selections{s}]),
        Mnull = mean(PSD.(dtype{d}).null.(selections{s}),2);
%         Enull = std(PSD.(dtype{d}).null.(selections{s}),0,2);
%         plot(freqs{s}, Mnull, 'k', 'linewidth', 2)
%         hold on,
%         plot(freqs{s}, Mnull+Enull, 'color', [0.3 0.3 0.3])
%         plot(freqs{s}, Mnull-Enull, 'color', [0.3 0.3 0.3])
        plot(PSD.(dtype{d}).freqs{s}, PSD.(dtype{d}).(selections{s})./Mnull, 'color', colors{s});
        set(gca, 'yscale', 'log', 'xscale', 'log');
        set(gca, 'ylim', [0.5 250], 'xlim', PSD.(dtype{d}).freqs{s}([1 end]));
        xlabel('frequency (Hz)'); ylabel('power ratio vs. trial shuffle')
    end
end
end

function [PSD, freqs] = makePSD(e, stimuli)
eNoise = e;
for stim = 1:8
    sel = find(stimuli==stim);
    meanE(:,stim) = smooth(mean(eNoise(:,sel),2), 0.03, 'lowess');
    eNoise(:,sel) = eNoise(:,sel) - meanE(:,stim);
end

for m =size(eNoise,2):-1:1
   [PSD(:,m), freqs] = pwelch(eNoise(:,m), 500,[],[], 1016); %PSD estimate
end
PSD = mean(PSD,2);
end

function summary = figure_PCA_perStim(dataset, resp, stimPeriod, prePeriod, postPeriod, doSub) %angles between first PCA components for ensembles in each stimulus direction, and spontaneous 
%for each segment, subtract the stimulus-associated activity
    respSub = resp;
    for stim = 1:8
        sel = find(dataset.stimulus.stim==stim);
        mStim(:,:,stim) = nanmean(respSub(:,:,sel),3);
        if doSub
            respSub(:,:,sel) = respSub(:,:,sel) - mStim(:,:,stim);
        end
    end
    %after mean subtraction, 
    %are the axes of variation for the different stimulus directions colinear? %(yes)
    nSVD = 10;
    w = nan(size(resp,1), nSVD, 9);
    s = nan(nSVD,9);
    for stim = 1:8
        sel = find(dataset.stimulus.stim==stim);
        D = reshape(respSub(:, stimPeriod,sel), size(respSub,1),[]); 
        ME = reshape(dataset.motion.error(stimPeriod,sel),1,[]);
        T = reshape(repmat(stimPeriod', 1, length(sel)), 1,[]);
        for shiftComp = 1:3  %3 motion axes obtained earlier by PCA of the 4-axis shifts
            M{stim}(:,shiftComp) = reshape(dataset.motion.shifts(shiftComp, stimPeriod,sel),[],1);
        end
        sel2 = any(isnan(D),1);
        D(:,sel2) = []; T(:,sel2) = []; ME(:,sel2) = []; M{stim}(sel2,:) = [];
        [U,S,V] = svd(D,'econ');
        w(:,:, stim) = U(:,1:nSVD)./sign(mean(U(:,1:nSVD),1));
        s(:,stim) = diag(S(1:nSVD,1:nSVD))./sqrt(size(D,2));
        h{stim} = V(:,1:nSVD)./sign(mean(U(:,1:nSVD),1));
        me{stim} = ME'; %motion
        t{stim} = T';
        sumD{stim} = sum(D,1)';
        
    end
    D = reshape(respSub(:, [prePeriod postPeriod],:), size(respSub,1),[]);
    ME = reshape(dataset.motion.error([prePeriod postPeriod],:),1,[]);
    sel2 = any(isnan(D),1);
    D(:,sel2) = []; ME(:,sel2) = [];
    for shiftComp = 1:3  %3 motion axes obtained earlier by PCA of the 4-axis shifts
            M{9}(:,shiftComp) = reshape(dataset.motion.shifts(shiftComp, [prePeriod postPeriod],:),[],1);
    end
    M{9}(sel2,:) = [];
    [U,S,V] = svd(D,'econ');
    w(:,:,9) = U(:,1:nSVD);
    s(:,9) = diag(S(1:nSVD,1:nSVD))./sqrt(size(D,2));
    h{9} = V(:,1:nSVD)./sign(mean(U(:,1:nSVD),1));
    me{9} = ME';
    sumD{9} = sum(D,1)';
    
    %cluster the components by similarity
    X = reshape(permute(w, [1 3 2]), size(w,1),[]);
    Xc = corr(X).^2; 
    
    %reorder components by their similarity
    for iter = 1:3
    for stimIX = 1:9
        notTaken = 1:nSVD;
        perm = nan(1,nSVD);
        for spontComp = 1:nSVD
            [~,maxIX] = max(mean(Xc(9*(notTaken-1)+stimIX, ((spontComp-1)*9)+[1:stimIX-1 stimIX+1:9]),2));
            perm(spontComp) = notTaken(maxIX); notTaken(maxIX) = [];
        end
        tt = stimIX:9:9*nSVD;
        Xc(tt, :) = Xc(tt(perm),:); Xc(:,tt) = Xc(:,tt(perm)); 
        X(:, tt) =  X(:, tt(perm));
        w(:,:,stim) = w(:,perm,stim);
        s(:,stim) = s(perm,stim);
        h{stim} = h{stim}(:,perm);
    end
    end
    

    Xd = sqrt(2*(1-sqrt(Xc)));
    labels = repmat(1:9, 1, length(Xc)/9); perplexity = 5;
    Y = []; cost = inf;
    figure('name', 't_sne visualization')
    for repeat = 1:5
        [Ytmp,Ctmp] = tsne_d(Xd, labels, 2, perplexity); %TSNE
        if Ctmp<cost
            cost = Ctmp;
            Y = Ytmp;
        end
    end
    sizes = (reshape(s', 1,[]).*100./max(s(:)));
    colors = repmat([hsv(8); 0 0 0], length(sizes)/9,1);
    figure('name', ['Best t-sne; cost: ' num2str(cost)]), scatter(Y(:,1),Y(:,2),sizes, colors, 'filled');
    figure, imagesc(Xc)
    
    %correlation of SVD components to various variables: total intensity,
    %sample movement, grating phase/time within stimulus?
    for cNum = nSVD:-1:1
        for stim = 9:-1:1
            corrToSumD(cNum,stim) = corr(sumD{stim}, h{stim}(:,cNum));
            for shiftComp = 1:3
                corrToMotion(cNum,stim,shiftComp) = corr(M{stim}(:,shiftComp), h{stim}(:,cNum));
            end
            corrToMotion(cNum,stim,4) = corr(sqrt(sum(M{stim}.^2,2)), h{stim}(:,cNum));
            if stim<9
                corrToTime(cNum,stim) = corr(t{stim}, h{stim}(:,cNum));
            end
        end
    end
    
    %make example figure
    figure, imagesc(Xc(1:8*9, 1:8*9));
    
%     figmodes = [1 2 3 5 10];
%     pos = [0 1.2 2.4 4 5.6];
%     figure('name', 'spatial correlations of modes')
%     for mode1 = 1:length(figmodes)
%         for mode2 = 1:length(figmodes)
%             hold on,
%             imagesc(pos(mode1) +[0 1], -(pos(mode2) +[0 1]), Xc((figmodes(mode1)-1)*9 + (1:9), (figmodes(mode2)-1)*9 + (1:9)))
%         end
%     end
%     axis image

figure,
figmodes = 1:8;
plot(0*reshape([corrToMotion(figmodes,:,1) nan(length(figmodes),1)]',1,[]), 'k', 'linewidth', 2)
hold on,
plot(reshape([corrToMotion(figmodes,:,1).^2 nan(length(figmodes),1)]',1,[]), 'linestyle', 'none', 'marker', 'o', 'markerfacecolor', [0.7 0.7 0.7], 'markeredgecolor', 'none', 'markersize', 10)
ylim([-0.001 0.1])
figure,
figmodes = 1:8;
plot(0*reshape([corrToMotion(figmodes,:,1) nan(length(figmodes),1)]',1,[]), 'k', 'linewidth', 2)
hold on,
plot(reshape([corrToMotion(figmodes,:,1) nan(length(figmodes),1)]',1,[]), 'linestyle', 'none', 'marker', 'o', 'markerfacecolor', [0.7 0.7 0.7], 'markeredgecolor', 'none', 'markersize', 10)
ylim([-0.3 0.3])
    
    
%     figure('name', 'colored by corr^2 to global photons');
%     cmp = squeeze(hsv2rgb(zeros(1,101), linspace(0,1,101), linspace(0,1,101)));
%     tmp =  reshape((corrToSumD').^2,1,[]); colors = cmp(round(100*tmp)+1,:);
%     scatter(Y(:,1),Y(:,2),sizes, colors, 'filled'); colormap(cmp); colorbar
%     
%     figure('name', 'colored by corr^2 to time');
%     cmp = squeeze(hsv2rgb(0.3*ones(1,101), linspace(0,1,101), linspace(0,1,101)));
%     tmp =  reshape((corrToTime').^2,1,[]); colors = cmp(round(100*tmp)+1,:);
%     select = repmat([true(8,1) ; false], nSVD,1);
%     scatter(Y(select,1),Y(select,2),sizes(select), colors, 'filled'); colormap(cmp); colorbar
%     
%     for shiftComp = 1:4
%         figure('name', ['colored by corr^2 to motion axis ' int2str(shiftComp)]);
%         cmp = squeeze(hsv2rgb((0.45+0.15*shiftComp)*ones(1,101), linspace(0,1,101), linspace(0,1,101)));
%         tmp =  reshape((corrToMotion(:,:,shiftComp)'.^2),1,[]); colors = cmp(round(100*tmp)+1,:);
%         scatter(Y(:,1),Y(:,2),sizes, colors, 'filled'); colormap(cmp); colorbar
%     end
    
    summary.Xc = Xc;
    summary.varExplained = sizes;
    summary.corrToMotion = corrToMotion;
    summary.corrToTime = corrToTime;
    summary.tSne = Y;
    summary.corrToSumD = corrToSumD;
end


function summary = figure_PCA_global(dataset, resp, stimPeriod, prePeriod, postPeriod, doSub, segSel) %angles between first PCA components for ensembles in each stimulus direction, and spontaneous 
%for each segment, subtract the stimulus-associated activity

    doNMF = false;

    respSub = resp(:,50:end-50,:);
    for stim = 1:8
        sel = find(dataset.stimulus.stim==stim);
        mStim(:,:,stim) = nanmean(respSub(:,:,sel),3) - nanmean(nanmean(respSub(:,[prePeriod]-49,sel),3),2);
        if doSub
            respSub(:,:,sel) = respSub(:,:,sel) - mStim(:,:,stim);
        end
    end
    if ~doNMF
        respSub = respSub-mean(respSub,2);
    end
    
    SL = 9*ones(size(respSub,2), size(respSub,3));  %stimlabel
    for stim = 1:8
        sel = find(dataset.stimulus.stim==stim);
        SL(stimPeriod-49,sel) = stim; 
    end
    
    nSVD = 5;
    D = reshape(respSub, size(respSub,1),[]);
    SP = zeros(size(respSub,2), size(respSub,3)); SP(stimPeriod-49,:) = 1; SP = reshape(SP, 1,[]);
    %ME = reshape(dataset.motion.error,1,[]);
    SL = reshape(SL, 1,[]);
    for shiftComp = 1:3  %3 motion axes obtained earlier by PCA of the 4-axis shifts
        M(:,shiftComp) = reshape(dataset.motion.shifts(shiftComp, 50:end-50,:),[],1);
    end
    T = reshape(repmat(1:size(respSub,2), 1,1,size(respSub,3)), 1,[]);
    sel2 = any(isnan(D),1);
    D(:,sel2) = []; T(:,sel2) = []; SP(:,sel2) = []; SL(:,sel2) = []; M(sel2,:) = [];
    %M(:,sel2) = [];
    
%     if doNMF
%         %NMF
%         [W,H] = nnmf(D,nSVD);
%         w = W; h= H';
%         globalComp = W*H;
%     else
        %SVD
        [U,S,V] = svd(D,'econ');
        w= U(:,1:nSVD)./sign(mean(U(:,1:nSVD),1));
        s= diag(S(1:nSVD,1:nSVD))./sqrt(size(D,2));
        h = V(:,1:nSVD)./sign(mean(U(:,1:nSVD),1));
%     end
    %m{stim} = M; %motion
    t = T'; sp = SP'; sl = SL';
    sumD = sum(D,1)';

    %correlation of SVD components to various variables: total intensity,
    %sample movement, grating phase/time within stimulus?
    for cNum = 1:nSVD
            summary.corrToSumD(cNum) = corr(sumD, h(:,cNum)).^2;
            summary.corrToTime(cNum) = corr(t, h(:,cNum)).^2;
            summary.corrToStimPeriod(cNum) = corr(sp, h(:,cNum)).^2;
            for shiftComp = 1:3
                summary.corrToMotion(cNum,shiftComp) = corr(M(:,shiftComp), h(:,cNum)).^2;
            end
            summary.corrToMotion(cNum,4) = corr(sqrt(sum(M.^2,2)), h(:,cNum)).^2;
            for stim = 1:9
                summary.corrToStim(cNum,stim) = corr(sl==stim, h(:,cNum)).^2;
            end
    end
    [~,gNum] = max(summary.corrToSumD);
    globalComp = U(:,gNum)*s(gNum)*V(:,gNum)';
    
    figure, tsel = 1:4000;
    subplot(2,1,1), plot(h(tsel,gNum)); xlabel('global component');
    subplot(2,1,2), plot(M(tsel,1)); xlabel('Motion axis 1')
%     subplot(4,1,3), plot(M(tsel,2)); xlabel('Motion axis 2')
%     subplot(4,1,4), plot(M(tsel,3)); xlabel('Motion axis 3')
    
    %images of these components
%     for comp = 1:nSVD
%         colorSegs(dataset, segSel, w(:,comp));
%         set(gcf, 'name', ['SVD component ' int2str(comp)]);
%     end
    
    %compare mean amplitude of global component vs amplitude of
    %stimulus-modulation
    if doSub
        summary.normGlobalAll = sum(var(globalComp,[],2),1);  %var(reshape(sum(globalComp,1),1,[]));
        for stim = 1:9
            summary.normGlobalStim(stim) = sum(var(globalComp(:, sl==stim),[],2)); %var(reshape(sum(globalComp(:, sl==stim),1),1,[]));;
        end
        summary.normStimAll = squeeze(sum(var(mStim,[],2),1));
        %normStimMod = squeeze(sum(mean(mStim(:, stimPeriod-50,:).^2,2),1)); %var(squeeze(sum(mStim,1))); 
    end
    
    %generate the global predictor
    %take the profile across segments, crop out nans, and use \ to identify the global predictor in time
    D2 = reshape(respSub, size(respSub,1),[]);
    sel3 =  ~any(isnan(D2),2); D2 = D2(sel3,:);
    wPred = w(sel3,:);
    globalPred= wPred\D2;
    globalPred = reshape(globalPred(gNum,:), size(respSub,2), size(respSub,3));
    
    summary.G = globalPred;
end

function summary = figure_NMF_global(dataset, resp, stimPeriod, prePeriod, postPeriod, segSel, dist) %angles between first PCA components for ensembles in each stimulus direction, and spontaneous 

% spontOnly = true; %Use only spontaneous data in defining global modes?
% if spontOnly
%     D =  resp(:,[prePeriod postPeriod],:);
% else
%     D = reshape(resp(:,50:end-50,:), size(resp,1),[]);
% end

spontOnly = true;

D = resp(:,50:end-50,:);
discardTrials =squeeze(max(sum(isnan(D),1),[],2))>3;
Dtmp = D(:,:,~discardTrials);
discardSegs = any(any(isnan(Dtmp),2),3); clear Dtmp;
D = D(~discardSegs,:,:);
summary.selectedSegs = ~discardSegs;
discardTrials =squeeze(max(sum(isnan(D),1),[],2))>0;
D = D(:,:,~discardTrials);
summary.selectedTrials = ~discardTrials;
summary.stim = dataset.stimulus.stim(~discardTrials);

% discardTrials =sum(isnan(D),1)>3; %we will discard some timepoints and some segments to get a data matrix with no nans
% Dtmp = D(:,~discardT);
% discardSegs = any(isnan(Dtmp),2); clear Dtmp;
% D = D(~discardSegs,:,:);
% summary.selectedSegs = ~discardSegs;
% discardT =sum(isnan(D),1)>0;
% D = D(:,~discardT);
rng(42); %fix the rng; this was done after the figures were made
nNMF = 8;
nShuffle=  9;
for repeat = 0:nShuffle
    if repeat
        for seg = 1:size(D,1)
            for stim = 1:8
                stimSel = find(summary(1).stim==stim);
                p = randperm(length(stimSel));
                D(seg,:,stimSel) =  D(seg,:,stimSel(p));
            end
        end
    end

    if spontOnly
        D2 = reshape(D(:,[prePeriod postPeriod]-49,:), size(D,1),[]);
    else
        D2 = reshape(D, size(D,1),[]);
    end
    opt = statset('MaxIter',5,'Display','final');
    [W0,H0] = nnmf(D2,nNMF,'replicates',10,'options',opt,'algorithm','mult');
    [summary(repeat+1).W,summary(repeat+1).H] = nnmf(D2,nNMF,'w0',W0,'h0',H0);

    if repeat<2
        %images of components
        sel = segSel;
        f = find(sel);
        sel(f(discardSegs)) = false;
%         for comp = 1:nNMF
%             colorSegs(dataset, sel, summary(repeat+1).W(:,comp));
%             set(gcf, 'name', ['SVD component ' int2str(comp) '. repeat: ' int2str(repeat)]);
%         end
        colorSegs(dataset, sel, summary(repeat+1).W(:,1:3));
        set(gcf, 'name', ['top three NMF components' '. Repeat: ' int2str(repeat)]);
    end

%Measure of distance scale:
binedges = 0:5:100;
bincenters = binedges+binedges(2)/2;
dd = dist(summary(1).selectedSegs,summary(1).selectedSegs);
Wnorm = summary(repeat+1).W./sqrt(sum(summary(repeat+1).W.^2,2));
Wnorm(isinf(Wnorm)|isnan(Wnorm)) = 0;
for bin = length(binedges):-1:1
    if bin==length(binedges)
        select = dd>binedges(bin);
    else
        select = dd>binedges(bin) & dd<=binedges(bin+1);
    end
    [ii,jj] = find(select);
    if any(ii)
    for comp = 1:nNMF
        summary(repeat+1).DCVbinned(comp,bin) =  corr(summary(repeat+1).W(ii,comp), summary(repeat+1).W(jj,comp));
        summary(repeat+1).DCVbinnedNorm(comp,bin) =  corr(Wnorm(ii,comp), Wnorm(jj,comp));
    end
    else
        summary(repeat+1).DCVbinned(comp,bin) = nan;
        summary(repeat+1).DCVbinnedNorm(comp,bin) = nan;
    end
end
%figure('name', ['Repeat: ' int2str(repeat)]), plot(bincenters(2:end), summary(repeat+1).DCVbinned(:,(2:end))')
%figure('name', ['Repeat: ' int2str(repeat) ' NORMALIZED']), plot(bincenters(2:end), summary(repeat+1).DCVbinnedNorm(:,(2:end))')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%use regression to solve for each mode on full data (mean-subtracted?)
doSub = true;
D2 = D;
if doSub
    for stim = 1:8
        sel = find(summary(1).stim==stim);
        mStim(:,:,stim) = nanmean(D2(:,:,sel),3);
        D2(:,:,sel) = D2(:,:,sel) - mStim(:,:,stim);
    end
end

modes = summary(repeat+1).W\reshape(D2,size(D2,1),[]);
modes = reshape(modes,  size(modes,1), size(D2,2), size(D2,3));
summary(repeat+1).modes = modes;

normfactor = sum(summary(repeat+1).W(:,1:3).^2,1)'; normfactor = normfactor./max(normfactor);
trialset=  2;
figure, plot(0.01*[1:3] + reshape(modes(1:3,:,8*(trialset-1)+(1:8)),3,[])'.*normfactor' )
stimstart = 950 + size(modes,2)*(0:7);
stimend = 3150 + size(modes,2)*(0:7);
stimX = [stimstart ; stimend ; nan*stimstart];
hold on, plot(stimX(:), 0.04*ones(size(stimX(:))), 'linewidth', 2)


for trial = size(modes,3):-1:1
    trial
    for mode1 = 3:-1:1
        for mode2 = 3:-1:1
            xcm(:,mode1,mode2,trial) = xcorr(modes(mode1,:,trial), modes(mode2,:,trial), 1000);
        end
    end
end
            


%tuning of each mode; calculate Z-scored orientation index
%relative power in each mode
summary(repeat+1).TC = []; summary(repeat+1).TC_E = [];

TCdata = var(modes(:,stimPeriod,:),[],2) - var(modes(:,prePeriod,:),[],2);
for stim = 1:8
    sel = summary(1).stim==stim;
    summary(repeat+1).TC(:,stim) = mean(TCdata(:,:,sel),3); 
    summary(repeat+1).TC_E(stim,:) = std(TCdata(:,:,sel),0,3)./sqrt(sum(sel));
end
for ORIshuff = 1000:-1:1
    TCdata = TCdata(:,:,randperm(size(TCdata,3)));
    for stim = 1:8
        sel = summary(1).stim==stim;
        summary(repeat+1).TCshuff(:,stim, ORIshuff) = mean(TCdata(:,:,sel),3);
        summary(repeat+1).TCshuff(:,stim, ORIshuff) = std(TCdata(:,:,sel),0,3)./sqrt(sum(sel));
    end
end
    

%frequency spectrum for each mode
% for mode = nNMF:-1:1
%     PSD(:,mode)  = makePSD(squeeze(modes(mode,:,:)), dataset.stimulus.stim);
% end
end
end

function summary = figure_ensembles(dataset, resp, stimPeriod, prePeriod, postPeriod)
respCut = resp(:,50:end-50,:);
respCat = reshape(respCut, size(respCut,1), []);

%we define a segment as active if its brightness is more than Q standard deviations of
%the trimmed (<80th percentile) distribution above its mean.
Q=3.5; %number of standard deviations above mean to call an event;
ampTMP = respCat;
ampTMP(ampTMP>prctile(ampTMP, 80,2)) = nan;
thresh = max(1.5, nanmean(ampTMP,2)+Q*nanstd(ampTMP,0,2));
[e, c, E] = computeEnsembles(respCut,thresh); %ensemble sizes and correlation matrix
binsE = 1:size(resp,1);
binsC = linspace(-1,1,200);
histE = histc(e(:), binsE); %histogram of ensemble size
histEspont = histc(reshape(e([prePeriod postPeriod]-50,:), 1, []), binsE);
histEevoked = histc(reshape(e([stimPeriod]-50,:), 1, []), binsE);
histC = histc(c(tril(true(size(c)),-1)), binsC); %histogram of pairwise correlations

% %frequency spectrum of global events
eNoise = e;
for stim = 1:8
    sel = find(dataset.stimulus.stim==stim);
    meanE(:,stim) = smooth(mean(eNoise(:,sel),2), 0.03, 'lowess');
    eNoise(:,sel) = eNoise(:,sel) - meanE(:,stim);
end

%shuffle trials per stimulus
doShuffle = 80;
eShuffS = nan([size(e) doShuffle]);
cShuffS = nan([size(c) doShuffle]);
for i = 1:doShuffle
    disp(['trial shuffle iteration ' int2str(i)])
    respTMP = respCut;
    for seg = 1:size(resp,1)
        for stim = 1:8
            sel = find(dataset.stimulus.stim==stim);
            respTMP(seg,:,sel) = respTMP(seg,:,sel(randperm(length(sel))));
        end
    end
    [eShuffS(:,:,i), cShuffS(:,:,i)] = computeEnsembles(respTMP,thresh);
end

histEshuff = histc(eShuffS(:), binsE);
histEshuff = histEshuff.* (sum(histE)./sum(histEshuff));
figure('name', 'ensemble size,ALL')
scatter((binsE(1:end-1)+binsE(2:end))/2,histE(1:end-1), 'markeredgecolor', [0 0 0])
hold on, plot((binsE(1:end-1)+binsE(2:end))/2, histEshuff(1:end-1), 'linewidth', 2, 'color', [0.5 0.5 0.5])
xlim([0 0.5*size(resp,1)]); ylim([1e-1 1e6]);
set(gca, 'xtick', [linspace(0, 0.8*size(resp,1), 9)], 'xticklabel', [0:0.1:0.8], 'YScale', 'log', 'tickdir', 'out')
xlabel('Ensemble Size (Fraction of Segments Recorded)')
ylabel('Frequency')
legend({'Observed', 'Trials shuffled within each condition'})

histEshuffSpont = histc(reshape(eShuffS([prePeriod postPeriod]-50, :,:),1,[]), binsE);
histEshuffSpont = histEshuffSpont.* (sum(histEspont)./sum(histEshuffSpont));
figure('name', 'ensemble size, SPONTANEOUS')
scatter((binsE(1:end-1)+binsE(2:end))/2,histEspont(1:end-1), 'markeredgecolor', [0 0 0])
hold on, plot((binsE(1:end-1)+binsE(2:end))/2, histEshuffSpont(1:end-1), 'linewidth', 2, 'color', [0.5 0.5 0.5])
xlim([0 0.5*size(resp,1)]); ylim([1e-1 1e6]);
set(gca, 'xtick', [linspace(0, 0.8*size(resp,1), 9)], 'xticklabel', [0:0.1:0.8], 'YScale', 'log', 'tickdir', 'out')
xlabel('Ensemble Size (Fraction of Segments Recorded)')
ylabel('Frequency')
legend({'Observed', 'Trials shuffled within each condition'})

histEshuffEvoked = histc(reshape(eShuffS([stimPeriod]-50, :,:),1,[]), binsE);
histEshuffEvoked = histEshuffEvoked.* (sum(histEevoked)./sum(histEshuffEvoked));
figure('name', 'ensemble size, EVOKED')
scatter((binsE(1:end-1)+binsE(2:end))/2,histEevoked(1:end-1), 'markeredgecolor', [0 0 0])
hold on, plot((binsE(1:end-1)+binsE(2:end))/2, histEshuffEvoked(1:end-1), 'linewidth', 2, 'color', [0.5 0.5 0.5])
xlim([0 0.5*size(resp,1)]); ylim([1e-1 1e6]);
set(gca, 'xtick', [linspace(0, 0.8*size(resp,1), 9)], 'xticklabel', [0:0.1:0.8], 'YScale', 'log', 'tickdir', 'out')
xlabel('Ensemble Size (Fraction of Segments Recorded)')
ylabel('Frequency')
legend({'Observed', 'Trials shuffled within each condition'})

%global tuning of dFF and #active segments
colors = hsv(8);
gDFF= squeeze(nansum(dataset.dPhotons, 1));
for stim = 1:8
    sel = find(dataset.stimulus.stim==stim);
    tmp = mean(gDFF(stimPeriod,sel),1) - mean(gDFF(prePeriod,sel),1);
    gTCdff_mean(stim) = mean(tmp);  gTCdff_err(stim) = std(tmp)./sqrt(sum(~isnan(tmp)));
    
    tmp  = (mean(e(stimPeriod-50,sel),1) - mean(e(prePeriod-50,sel),1))./size(resp,1);
    gTCe_mean(stim) = mean(tmp);  gTCe_err(stim) = std(tmp)./sqrt(sum(~isnan(tmp)));
end
figure, 
subplot(2,1,1), errorbar(gTCdff_mean, gTCdff_err, 'color', 'k', 'linewidth', 2, 'linestyle', ':');
xlabel('stimulus'); ylabel('dPhotons/frame'); set(gca, 'xlim', [0 9]);
subplot(2,1,2), errorbar(gTCe_mean, gTCe_err, 'color', 'b', 'linewidth', 2, 'linestyle', ':');
xlabel('stimulus'); ylabel('#active segments/frame'); set(gca, 'xlim', [0 9]);
for stim = 1:8
    subplot(2,1,1)
    hold on, scatter(stim, gTCdff_mean(stim), 'sizedata', 120, 'markerfacecolor', colors(stim,:), 'markeredgecolor', 'k', 'linewidth', 2)
    subplot(2,1,2)
    hold on, scatter(stim, gTCe_mean(stim), 'sizedata', 120, 'markerfacecolor', colors(stim,:), 'markeredgecolor', 'k', 'linewidth', 2)
end

%frequency of events exceeding the 99% confidence limit of the null, per stimulus
prctCutoff= 99;
for stim = 1:8
      sel = find(dataset.stimulus.stim==stim);
      TCthresh(stim) = prctile(reshape(eShuffS(stimPeriod,sel,:), 1,[]), prctCutoff);
      TC(stim) = mean(reshape(e(stimPeriod,sel)>TCthresh(stim),1,[]));
      N(stim) = length(stimPeriod)*length(sel);
end
TCthresh(9) = prctile(reshape(eShuffS([prePeriod postPeriod]-50,:,:), 1,[]),  prctCutoff);
TC(9) = mean(reshape(e([prePeriod postPeriod]-50,:)>TCthresh(9),1,[]));
N(9) = length([prePeriod postPeriod])*size(e,2);

summary.TC = TC;
summary.TCN = N;
summary.e = e;
summary.eShuff = eShuffS;
summary.evoked.hist = histEevoked;
summary.evoked.histShuff = histEshuffEvoked;
summary.spont.hist = histEspont;
summary.spont.histShuff = histEshuffSpont;
summary.bins = binsE;
end

function summary = figure_crossCorrelations(resp, d, label, dX,dY,stimuli)
if nargin<3
    label = '';
end
    %correlate spikes
    binedges = 0:5:100;
    maxlags = 700;
    summary.XC = zeros(maxlags+1, size(resp,1).^2);
    summary.XCx = zeros(maxlags+1, length(binedges), max(stimuli));
    summary.XCy = zeros(maxlags+1, length(binedges), max(stimuli));
    summary.XCn = zeros(length(binedges), max(stimuli));
    N = zeros(1, size(summary.XC,2));
    for m = 1:size(resp,3)
        m
        %Total crosscorrelation
        tmp = xcorr((resp(:,:,m)-mean(resp(:,:,m),2))', maxlags, 'none');
        tmp = tmp(maxlags+1:end,:);
        sel = ~any(isnan(tmp),1);
        summary.XC(:,sel) = summary.XC(:,sel) + tmp(:,sel);
        N(sel) = N(sel)+1;
        
        %separated by stimuli
        stimID = stimuli(m);
        for bin = 1:length(binedges)
            if bin==length(binedges)
                selectX = dX>binedges(bin);
                selectY = dY>binedges(bin);
            else
                selectX = dX>binedges(bin) & dX<=binedges(bin+1);
                selectY = dY>binedges(bin) & dY<=binedges(bin+1);
            end
            summary.XCn(bin,stimID) = summary.XCn(bin,stimID) + sum(~isnan(tmp(1, selectX)));
            corr0 = reshape(tmp(1,:), size(resp,1), size(resp,1));
            normfactor = sqrt(diag(corr0).*diag(corr0)');
            
            summary.XCx(:,bin,stimID) = summary.XCx(:,bin,stimID) + nansum(tmp(:,selectX),2)./nansum(normfactor(selectX));       %delays,bin, stimulus
            summary.XCy(:,bin,stimID) = summary.XCy(:,bin,stimID) + nansum(tmp(:,selectY),2)./nansum(normfactor(selectY));       %delays,bin, stimulus
        end
    end
    summary.XC = summary.XC./N;
    [~, summary.maxCorrDelayX] = max(summary.XCx,[],1); summary.maxCorrDelayX = squeeze(summary.maxCorrDelayX);
    [~, summary.maxCorrDelayY] = max(summary.XCy,[],1); summary.maxCorrDelayY = squeeze(summary.maxCorrDelayY);
    
    meanCorr0 = reshape(summary.XC(1,:), size(resp,1), size(resp,1));
    normfactor = sqrt(diag(meanCorr0).*diag(meanCorr0)');
    summary.maxCorrDelay = nan(1,length(binedges));
    summary.XCbinned = nan(maxlags+1, length(binedges));
    summary.binedges = binedges;
    for bin = 1:length(binedges)
        if bin==length(binedges)
            select = d>binedges(bin);
        else
            select = d>binedges(bin) & d<=binedges(bin+1);
        end
        tmp = summary.XC(:, select);
        summary.XCbinned(:,bin) =  nanmean(tmp,2)./nanmean(normfactor(select));
        [~,maxix] = max(summary.XCbinned(:,bin));
        summary.maxCorrDelay(bin) = maxix;
    end
    figure('name', label), imagesc([0 size(summary.XCbinned,1)/1.016], [(binedges(2)+binedges(3))/2 (binedges(end)+binedges(2)/2)], summary.XCbinned(:,2:end)')
    hold on, scatter(summary.maxCorrDelay, binedges+(binedges(2)/2), 'marker', '*', 'markeredgecolor', 'r');
    set(gca, 'ytick', [binedges(1:end) binedges(end)+binedges(2)], 'tickdir', 'out')
    xlabel('delay (ms)'); ylabel('distance (um)')
end

function [d,dX,dY] = find_distances(dataset,segSel)
nSeedSel = sum(segSel);
d = nan(nSeedSel);
dims = size(dataset.refIM.IM);
segs = dataset.refIM.seg(:,segSel);
segs = reshape(full(segs),[dims,nSeedSel]);
S = squeeze(any(segs,3));
S = reshape(S,[dims(1)*dims(2),nSeedSel]);

%% find segment centers
centers = nan(nSeedSel,2);

for seed = 1:nSeedSel
    [cdummyX, cdummyY] = ind2sub(dims(1:2),find(S(:,seed)));
    centers(seed,:) = [mean(cdummyX),mean(cdummyY)];
end
centers = dataset.refIM.segopts.XYscale*centers; %in um

disp('Calculating segment distances ...')
dX = abs(centers(:,1) - centers(:,1)');
dY = abs(centers(:,2) - centers(:,2)');
d = sqrt(dX.^2 + dY.^2); d(logical(eye(size(d)))) = nan;

disp('done')
if min(d(:))< 0.001
    disp('There are fused into each other segments!');
    keyboard;
end
unitDistance = 5; %um
nBins = 10;
Binned = false(nSeedSel,nSeedSel,nBins);
for bin = 1:nBins
    if bin ~= nBins
        Binned(:,:,bin) = d >= (bin-1)*unitDistance & d < bin*unitDistance;
    else
        Binned(:,:,bin) = d >= (bin-1)*unitDistance;
    end
end
end

function [ensembleSize, C, E] = computeEnsembles(resp,thresh)
    respCat = reshape(resp, size(resp,1), []);
    C = corr(respCat');
    E = resp>thresh;
    ensembleSize = squeeze(sum(E,1));
end

function [H,S,V] = circular_mean(rad,w)
A = sum(exp(1i*rad).*w,2);
H = 1/2 + angle(A)/(2*pi);

S = abs(A)./sum(abs(w),2);
V = sqrt(sum(w.^2,2));
end

function h2 = HPfilter(h)
    %h2 =imtophat(max(h,0), ones(1,50,1));
    
    dim = 2;
    cutOffFreqHz = 4;
    cutoff = cutOffFreqHz/1016;
    cutoffIX = round(cutoff*size(h,dim));
    H = fft(h,[], dim);
    H(:,[1:cutoffIX end-cutoffIX+1:end],:) = 0;
    h2 = real(ifft(H, [],dim));
end

function colorSegs(dataset, segSel, colors)
    S = dataset.refIM.seg(:, segSel);
    ncolors= size(colors,2);
    if ncolors==1
        IM = reshape(full(sum(S.*colors',2)), size(dataset.refIM.IM));
    else
        IM = zeros([size(dataset.refIM.IM) 3]);
        for c = 1:ncolors
            IMtmp =reshape(full(sum(S.*colors(:,c)',2)), size(dataset.refIM.IM));
            IMtmp = IMtmp./prctile(IMtmp(:),99.995);
            IM(:,:,:,c) = IMtmp;
        end
    end
    if size(IM,4)>1
        figure('name', 'Avg intensity projection'),
        IMtmp = squeeze(mean(flipdim(permute(IM, [2 1 3 4]),1),3)); IMtmp = IMtmp./prctile(IMtmp(:),99.995);
        subplot(1,2,1), imshow(IMtmp);
        subplot(1,2,2), imshow(mean(dataset.refIM.IM,3),[]);
    end
    
    figure('name', [int2str(sum(segSel)) 'selected segments']), imshow3D(flipdim(permute(IM, [2 1 3 4]),1));
end