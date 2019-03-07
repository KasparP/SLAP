function opts = defaultOpts_InVivoGluSnFr

opts.tau = 25;                  %'Decay time constant';
opts.dampar = 1.5;              %'Residuals less than dampar S.D.s will have less effect on updates; typically 0-2';
opts.numIters = 50;            %'Maximum number of iterations';
opts.minNumIters = 50;          %'Minimum number of iterations';
opts.darkNoise = 0.05;          %'Dark noise, in photons, added to all expected photon rates';
opts.stimArtifact = 0;          %'set to length of artifact to compensate for a widefield stimulus artifact';
opts.clipGradients = true;      %'Clip multiplicative updates to 0.1<X<10 per iteration';
opts.offFocusSpikes = true;     %Global spikes outside the focal plane; set to true if out-of-focus background can produce large coherent signals
opts.baselineMult = 1;          %Multiplier for the calculated baseline, typically 0.5-1; set to 1 if alignment is good, lower if alignment is poor
opts.correctBleaching = true;   %fit a rank-1 temporally varying baseline to capture global bleaching
opts.threshold_spikes = true;   %set small spikes to 0
opts.clipDFF = 3;               %maximum DFF, set to Inf to not clip
end