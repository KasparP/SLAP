function converged = check_convergence(sys)
dff_converged = false;
log_likelihood_coverged = false;
        if sys.opts.convergence_flag
            switch sys.opts.baseline_model
                case 'multiplicative'
            SLAPMi_messages('objective',sys)
            if log10(sys.control_params.dff_diff(sys.control_params.current_iter)) < log10(sys.opts.convergence_threshold)-3
                SLAPMi_messages('dff break',sys);
                dff_converged = true;
            end
                case 'additive_bRef'
                case 'additive_constant_baseline'
                otherwise
                    SLAPMi_messages('baseline_model',sys);
            end
            rl = 1;
            sys.window_length = 5;
            if sys.control_params.current_iter > sys.opts.min_num_iters
                indices = sys.control_params.current_iter-sys.window_length:sys.control_params.current_iter-1;
                base = interp1(indices,sys.control_params.GoF(indices),sys.opts.NumIter_dynamic,'linear','extrap');
                A = sys.control_params.GoF - base;
                rl = abs(1-A(indices)./A(indices+1));
            end
            ll_conv = log10(rl)< log10(sys.opts.convergence_threshold)-2;
            if sum(ll_conv) == sys.window_length
                SLAPMi_messages('log-likelihood break',sys);
                log_likelihood_coverged = true;
            end
        end
        
%         converged = (dff_converged || log_likelihood_coverged) & ~sys.solver_params.regularizatio_flag ;
        converged = (dff_converged || log_likelihood_coverged) ;
end

