function sys = initialize_solver (sys,opts)

[sys.opts.n, sys.opts.T] = size(sys.input.y);

sys.opts.p = size(sys.output.PS,2);

%USER OPTIONS
sys.opts.tau = opts.tau;
sys.opts.theta = [1 -exp(-1/sys.opts.tau)];
sys.opts.NumIter_dynamic = opts.numIters;
sys.opts.min_num_iters = opts.minNumIters;
sys.solver_params.dampar = opts.dampar;
sys.opts.clip_gradients = opts.clipGradients;
sys.opts.dark_noise = opts.darkNoise;
if ~isfield(sys.solver_params,'solver'); sys.solver_params.solver = 'dyn_RL';                  end
sys.opts.convergence_flag = true;
sys.opts.convergence_threshold = 1; %in percentage
if ~isfield(sys.solver_params,'verbose'); sys.solver_params.verbose = 100;                     end
sys.opts.ReconAlg = opts.ReconAlg;
sys.opts.offFocusSpike = opts.offFocusSpikes;

% Regularization
sys.solver_params.regularization_flag = false;
if ~isfield(sys.solver_params,'lambda');  sys.solver_params.lambda = 0;                        end
if ~isfield(sys.solver_params,'pen_norm'); sys.solver_params.pen_norm = 'vec_norm';            end
if ~isfield(sys.solver_params,'dff_pnorm'); sys.solver_params.dff_pnorm = 2;                   end
if ~isfield(sys.solver_params,'dff_qnorm'); sys.solver_params.dff_qnorm = 10;                  end
if ~isfield(sys.solver_params,'max_dff'); sys.solver_params.max_dff = 4;                       end

%BASELINE
sys.opts.baseline_model = 'additive_constant_baseline';
sys.opts.baselineMult = opts.baselineMult;

%Initialization
sys.opts.filter_delays = [];
sys.control_params.dff_diff = zeros(1,sys.opts.NumIter_dynamic);
sys.control_params.regularization = zeros(1,sys.opts.NumIter_dynamic);
sys.control_params.penalty = zeros(1,sys.opts.NumIter_dynamic);
sys.control_params.Likelihood = zeros(1,sys.opts.NumIter_dynamic);
sys.opts.momentum = 1;
sys.opts.upper_spike_bounds = sys.solver_params.max_dff*filter(sys.opts.theta,1,ones(sys.opts.p,sys.opts.T),sys.opts.filter_delays,2);
sys = initialize_spikes(sys);
% sys.output.spikes(end,:) = 0;

if sys.opts.correctBleaching
    y  = max(imgaussfilt(sys.input.y,[eps 20], 'padding', 'symmetric'), imgaussfilt(sys.input.y,[0.5 20], 'padding', 'symmetric'));
    sys.output.additive_baseline = max(min(y,[],2), sys.opts.dark_noise);
    btemp = sys.output.additive_baseline'*sys.input.y/(sys.output.additive_baseline'*sys.output.additive_baseline);
    btemp = smooth(cummin(imgaussfilt(btemp,[eps,11],'padding','symmetric')),500)';
    sys.output.baseline_temp = fit_spline(btemp);
    y2b = sum(sys.input.y)./sum(sys.output.additive_baseline*sys.output.baseline_temp);
    alpha = prctile(y2b(y2b>0),15);
    sys.output.baseline_temp = opts.baselineMult.*alpha*sys.output.baseline_temp;
else
    y  = max(imgaussfilt(sys.input.y,[eps 100], 'padding', 'symmetric'), imgaussfilt(sys.input.y,[0.5 100], 'padding', 'symmetric'));
    sys.output.additive_baseline = max(opts.baselineMult.*min(y,[],2), sys.opts.dark_noise);
end
sys = calc_rates(sys);

SLAPMi_messages('Time Constant',sys)
SLAPMi_messages('iterations',sys);
SLAPMi_messages('mdff',sys);
SLAPMi_messages('regularization',sys)
SLAPMi_messages('ref_image',sys);
end

function yOut = yRef(yIn)
    %get an estimate of the F0 line intensities
    yOut = imgaussfilt(yIn,[eps 50], 'padding', 'symmetric');
    yIn(isnan(yIn)) = nanmin(yOut(:));
    yOut = imgaussfilt(yIn,[eps 50], 'padding', 'symmetric');
    yOut = min(yOut,[],2);
end
