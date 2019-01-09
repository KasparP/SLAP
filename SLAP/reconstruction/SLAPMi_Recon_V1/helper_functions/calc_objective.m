function sys = calc_objective( sys )
%%% This function calculates the objective function
% at the current estimates
switch sys.solver_params.regularization_flag
    case true
        switch sys.solver_params.pen_norm
            case 'p_norm'
            sys.penalty(sys.control_params.current_iter) = sum(sum(sys.dff.^sys.solver_params.dff_pnorm,2).^(sys.solver_params.dff_qnorm/sys.solver_params.dff_pnorm))^(1/sys.solver_params.dff_qnorm);
            case 'vec_norm'
                sys.penalty(sys.control_params.current_iter) = sum(sys.dff(:).^sys.solver_params.dff_qnorm);
            otherwise
                SLAPMi_messages('regularization',sys);
        end
    case false
        sys.penalty(sys.control_params.current_iter) = 0;
end
switch sys.opts.ReconAlg
    case 'RL'
        L = sys.control_params.rate - sys.input.y.*log(sys.control_params.rate);
    case 'nnls'
        L = norm(sys.input.y-sys.control_params.rate,'fro')^2;
    otherwise
        SLAPMi_messages('ReconAlg',sys)
end
        sys.control_params.Likelihood(sys.control_params.current_iter) = sum(L(:));
        sys.control_params.objective(sys.control_params.current_iter) = sys.control_params.Likelihood(sys.control_params.current_iter) +...
            sys.solver_params.lambda*sys.control_params.penalty(sys.control_params.current_iter);
%   Goodness of fit 
        sys.control_params.GoF(sys.control_params.current_iter) = sys.control_params.Likelihood(sys.control_params.current_iter);
        
end

