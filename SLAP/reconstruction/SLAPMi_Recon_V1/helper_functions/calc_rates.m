function sys = calc_rates( sys )
% calc_rates(sys) calculates the estimated rates for the Poission lixel
% intensity model
if sys.opts.correctBleaching
sys.control_params.rate = sys.output.PS*sys.output.F+...
    sys.output.additive_baseline*sys.output.baseline_temp+sys.opts.dark_noise;
else
sys.control_params.rate = sys.output.PS*sys.output.F+sys.output.additive_baseline+sys.opts.dark_noise;
end

if sys.opts.removeStimArtifactFlag
    sys.control_params.rate = sys.control_params.rate + sys.output.stim_artifact;
end

% Calculate the variance of the estimates using Fisher information 
sys.output.FisherF = (sys.output.PS.^2)'*(1./sys.control_params.rate);
% approximately
sys.output.varF = 1./sys.output.FisherF;
end

