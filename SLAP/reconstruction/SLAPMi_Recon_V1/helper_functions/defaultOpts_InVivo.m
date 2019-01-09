function opts = defaultOpts_InVivo
opts.tau = 200;             %'Decay time constant';
opts.dampar = 0;          %'Residuals less than dampar S.D.s will have less effect on updates; typically 0-2';
opts.numIters = 200;        %'Maximum number of iterations';
opts.minNumIters = 30;      %'Minimum number of iterations';
opts.darkNoise = 0.01;      %'Dark noise, in photons, added to all expected rates';
opts.stimArtifact = 0;      %'set to length of artifact to compensate for a widefield stimulus artifact';
opts.clipGradients = true;  %'Clip multiplicative updates to 0.1<X<10 per iteration';
opts.offFocusSpikes = true;
opts.baselineMult = 0;
opts.perSeedBaseline = false;
end