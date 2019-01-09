function sys  = update_spikes( sys )
%%% This function calculates the multiplicative updates
%   for the spikes (called "f" here).
Gn     = max(-sys.Gp,0);
sys.Gp = max(sys.Gp,0);
R1 = ((sys.Ln+sys.solver_params.lambda*Gn)./(sys.Lp+sys.solver_params.lambda*sys.Gp));
R2 = ((sys.Ln2+sys.solver_params.lambda*Gn)./(sys.Lp2+sys.solver_params.lambda*sys.Gp));
if sys.solver_params.lambda > 0 %|| sys.control_params.current_iter<20 %we are regularizing, havev't tested the two models update in these conditions
    R = R1;
else %use the 'two models' update
    R = min(R1,R2);
    %             R = R1;
end
if sys.opts.clip_gradients
    R =  max(min(R,10 ),0.1);
end
sys.output.spikes = sys.output.spikes.* R;
sys.output.spikes(isnan(sys.output.spikes)) = 0;

sys.control_params.X_temp = sys.output.F;

% threshold spikes
if sys.opts.threshold_spikes
sys.output.spikes(sys.output.spikes<1e-10) = 0;
end

%   Update df/f or X
sys.output.F = filter(1,sys.opts.theta,sys.output.spikes,[],2);

%threshold DFF if necessary
if isfield(sys.opts, 'clipDFF') && sys.opts.clipDFF>0 && ~isinf(sys.opts.clipDFF)
    if any(sys.output.F(:)>sys.opts.clipDFF)
        F = sys.output.F;
        minF = min(F(:,50:end-50),[],2);
        F = (F-minF)./(1+minF);
        F = max(0, min(F, sys.opts.clipDFF)); % min(F, sys.opts.clipDFF)
        sys.output.F = (F.*(1+minF))+minF;
        spikes = filter(sys.opts.theta,1,sys.output.F,sys.opts.filter_delays,2);
        sys.output.spikes = max(0,spikes);
    end
end

sys.control_params.dff_diff(sys.control_params.current_iter) = ...
    norm(sys.output.F-sys.control_params.X_temp,'fro')/norm(sys.control_params.X_temp,'fro');



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%estimate per seed baseline when half of iterations are done:
%   Calculate the current value of the objective function
%   Can be used for sanity checks
%%  Update the current estimates of the rates ...
%   based on the estimated df/f
sys = calc_rates(sys);
sys = calc_objective(sys);
% if sys.control_params.current_iter == 20
%     sys.opts.ReconAlg = 'RL';
% end
end

