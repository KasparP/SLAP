function sweep_parameters(problemNumber)
    tau = [1,10];
    dampar = [0 1.5];
    offFocusSpikes = [true,false];
    baselineMult = [0.5 1];
    
    index = dec2bin(problemNumber,4);

    opts.tau = tau(str2double(index(1))+1);
    opts.dampar = dampar(str2double(index(2))+1);
    opts.numIters = 100;        %'Maximum number of iterations';
    opts.minNumIters = 50;      %'Minimum number of iterations';
    opts.darkNoise = 0.01;      %'Dark noise, in photons, added to all expected rates';
    opts.stimArtifact = 0;      %'set to length of artifact to compensate for a widefield stimulus artifact';
    opts.clipGradients = true;  %'Clip multiplicative updates to 0.1<X<10 per iteration';
    opts.offFocusSpikes = offFocusSpikes(str2double(index(3))+1);
    opts.baselineMult = baselineMult(str2double(index(4))+1);
    
    load fns
    opts.skipFrames = 2;
    opts.appendname = index;
    opts.plot = true;
    sys_recon = SLAPMi_solve(opts,{fns{[24,26,30,33]}});
end