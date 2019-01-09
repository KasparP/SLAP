function opts = defaultOpts_RedGreen

opts.tau = 50;             %'Decay time constant';
opts.dampar = 1.5;          %'Residuals less than dampar S.D.s will have less effect on updates; typically 0-2';
opts.numIters = 100;        %'Maximum number of iterations';
opts.minNumIters = 50;      %'Minimum number of iterations';
opts.darkNoise = 0.01;      %'Dark noise, in photons, added to all expected rates';
opts.stimArtifact = 0;      %'set to length of artifact to compensate for a widefield stimulus artifact';
opts.clipGradients = true;  %'Clip multiplicative updates to 0.1<X<10 per iteration';
opts.offFocusSpikes = false;
opts.baselineMult = 1;
opts.perSeedBaseline = false;
opts.ReconAlg = 'RL';
% opts.StimOnset = 500;
end