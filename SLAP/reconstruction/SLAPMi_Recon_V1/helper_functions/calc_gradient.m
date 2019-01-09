function sys = calc_gradient( sys )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function calculates the positive and negative components
% in the gradient of the objective function wrt the spikes

switch sys.opts.ReconAlg
    case 'RL'
        %% Step 1:
        %          Normalize the observations by the updated rates
        %          This only goes into calculation of the negative terms
        %          Please check the manuscript for details
        y_r = sys.input.y./sys.control_params.rate;
%         y_r = min(y_r,5);
%         y_r = dampFactor(sys.input.y,sys.control_params.rate);
        %          For the positive terms the corresponding step is
        %          equivalent to replacing the normalized measurements by 1
        %          Calculate the negative and positive terms wrt df/f
        if sys.solver_params.dampar == 0 || sys.control_params.current_iter < 5
            y_n = sys.output.PS'*y_r;
        else
            N = 10;
            U0 = 2/sys.solver_params.dampar^2*(sys.input.y.*log(y_r)-sys.input.y+sys.control_params.rate);
            U = min(U0,1);
            y_n = (sys.output.PS'* (1+ (y_r-1).* U.^(N-1).*(N - (N-1)*U)));
        end
        %% Step 2: Map the gradient terms from df/f to spikes
        % Positive terms in the gradient of log likelihood
        y_p = repmat(sum(sys.output.PS)',1,sys.opts.T);   
    case 'nnls'
        % no dampar is going to be in effect here
        y_n = sys.output.PS'*sys.input.y;
        y_p = sys.output.PS'*sys.output.PS*sys.output.F;
    otherwise
        SLAPMi_messages('ReconAlg',sys)
end

sys.Lp = fliplr(filter(1,sys.opts.theta,full(fliplr(y_p)),sys.opts.filter_delays,2))+eps;
sys.Lp2 = y_p;
% Negative terms in the log likelihood
sys.Ln  = fliplr(filter(1,sys.opts.theta,fliplr(y_n),sys.opts.filter_delays,2));
sys.Ln2 = y_n;

%% Step 3: Calculate the derivative of the penalty function
%       wrt the spikes. This term is (almost) always positive
%       and gets added to the postive terms from the log likelihood
if sys.solver_params.lambda == 0
    sys.Gp = 0;
else
    sys.Gp = calc_norm_gradient(sys);
    sys.control_params.penalty_flag = true;
end
end

function [M,V] = poisson_dev(lambda )
% Prob( X > C lambda) <= exp[lambda*(C-1-Cln(C))] Chernoff Bound
% For lambda = 1, C = 5 gives a probability of 0.01;
M = 2*exp(-lambda).* (lambda.^(floor(lambda)+1))./factorial(floor(lambda));
% M = 2*exp(-lambda).* (lambda.^(lambda+1)./gamma(lambda+1);
V = lambda - M.^2;
end

