function Gp = calc_norm_gradient( sys )

    switch sys.solver_params.pen_norm
        case 'p_norm'
            if ~isfield(sys.control_params,'penalty_flag')
                SLAPMi_messages('regularization',sys)
%                 sprintf('Penalizing the %d - %d norm of delta_f \n',sys.solver_params.dff_pnorm,sys.solver_params.dff_qnorm)
            end
            A = sys.dff.^(sys.solver_params.dff_pnorm-1);
            Gp =  fliplr(filter(1,sys.opts.theta,fliplr(A),[],2));
            B = (repmat(sum(sys.dff.^sys.solver_params.dff_pnorm,2),1,sys.opts.T)).^(sys.solver_params.dff_qnorm/sys.solver_params.dff_pnorm-1);
            B(isinf(B)) = 0;
            B(isnan(B)) = 0;
            C = sum(sum(sys.dff.^sys.solver_params.dff_pnorm,2)...
                .^(sys.solver_params.dff_qnorm/sys.solver_params.dff_pnorm))^(1/sys.solver_params.dff_qnorm-1);
            Gp = Gp.*B*C;

        case 'vec_norm'
            if ~isfield(sys.control_params,'penalty_flag')
                SLAPMi_messages('regularization',sys)
            end
%             A = sys.dff.^(sys.solver_params.dff_qnorm-1);
            A = sys.dff.^(sys.solver_params.dff_qnorm-1)*(sum(sys.dff(:).^sys.solver_params.dff_qnorm)^(1/sys.solver_params.dff_qnorm-1));
            Gp =  fliplr(filter(1,sys.opts.theta,fliplr(A),[],2));
        case 'infinity'
            if ~isfield(sys.control_params,'penalty_flag')
                disp('Only penalizing the maximum delta_f')
            end
%             [m,I] = max(sys.dff(:));
%             A = zeros(size(sys.dff));
%             A(I) = 1;
            I = sys.dff > sys.solver_params.max_dff;
            A = zeros(size(sys.dff));
            A(I) = sys.dff(I);
            Gp =  fliplr(filter(1,sys.opts.theta,fliplr(A),[],2));
        case 'nuclear'
            if ~isfield(sys.control_params,'penalty_flag')
                SLAPMi_messages('regularization',sys)
            end
            p = 1;
            gamma = 1e-4;
            W = (sys.output.F*sys.output.F' + gamma*eye(sys.opts.p))^(p/2-1);
            A = 2 * W * sys.output.F;
            Gp =  fliplr(filter(1,sys.opts.theta,fliplr(A),[],2));
    otherwise
        error('penalty norm not specified')
    end
    
% Gp(end,:) = 0;
% disp('No penalty on the background')
end

