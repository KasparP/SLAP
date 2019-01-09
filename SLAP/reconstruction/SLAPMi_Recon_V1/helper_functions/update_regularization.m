function  sys  = update_regularization( sys )
% A = sys.dff(:,end+floor(1/log(-sys.theta(2))));
% A = sys.dff(:,1:end-floor(sys.opts.tau/2));
A = sys.dff;
switch sys.solver_params.pen_norm
    case 'p_norm'
        m = sum(sum(A.^sys.solver_params.dff_pnorm,2)...
                .^(sys.solver_params.dff_qnorm/sys.solver_params.dff_pnorm))^(1/sys.solver_params.dff_qnorm);
    case 'vec_norm'
        m = max(A(:))/sys.solver_params.max_dff;
    otherwise
        SLAPMi_messages('penalty',sys);
end
        if m > 1
            if sys.solver_params.lambda == 0
                SLAPMi_messages('regularization',sys);
                sys.solver_params.lambda = 0.01;
%                 sys.solver_params.lambda = sys.solver_params.lambda*m;
%             else
%                 B = (sys.output.spikes.*sys.Ln./sys.opts.upper_spike_bounds-sys.Lp)./sys.Gp;
%                 sys.solver_params.lambda = max(B(:));
            end
        end
        
        sys.solver_params.lambda = sys.solver_params.lambda*m;
        
        
%         m=sqrt(m);
%         if (m-1)*(sys.opts.momentum-1) > 0 
%         sys.opts.momentum = sys.opts.momentum * m;
%         else
%         sys.opts.momentum = m;    
%         end
%         
%         sys.solver_params.lambda = sys.solver_params.lambda*sys.opts.momentum;
        
        sys.control_params.regularization(sys.control_params.current_iter) = sys.solver_params.lambda;
end