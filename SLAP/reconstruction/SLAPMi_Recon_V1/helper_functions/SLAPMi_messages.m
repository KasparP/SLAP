function SLAPMi_messages(mname,sys)
switch mname
    case 'negative_meas'
        disp('Measurement has negative values -- will be set to 0 !!')
        fprintf('Minimum of measurements is %d ...\n', min(sys.input.y(:)));
    case 'Time Constant Changed'
        if sys.control_params.current_iter == 1
         disp('Time constant changes from short to long over iterations')
        end
        
    case 'Time Constant'
        fprintf('Estimated time constant: %d ms \n',sys.opts.tau);
        
    case 'initialization'
        switch sys.initialization
        case 'all_ones'
        disp('Initilizing the solver with all one Delta F / F ...')
        case 'back_projection'  
        disp('Initilizing the solver with back projection solution ...')
        end
        
    case 'mdff'
        if sys.solver_params.regularization_flag
            switch sys.solver_params.max_dff
                case 20
                    sensor = 'GCAMP';
                case 4
                    sensor = 'yglusnfr';
                otherwise
                    sensor = 'unknown sensor';
            end
            fprintf('Maximum Delta F / F set to %d for %s ... \n',sys.solver_params.max_dff,sensor);
        end
        
    case 'iterations'
        fprintf('The solver is running for a maximum of %d iterations \n',sys.opts.NumIter_dynamic);
        
    case 'current iteration'
        if mod(sys.control_params.current_iter,sys.solver_params.verbose) == 0
        fprintf('Reached iteration %d ... \n',sys.control_params.current_iter);
        end
        
    case 'penalty'
        error('penalty type not specified')
        
    case 'stimArtifact'
        switch sys.opts.removeStimArtifactFlag
            case true
                disp('Stimulus artifact for two-spot uncaging dataset is removed ... ')
                fprintf('Stimuli are %d milliseconds apart \n',sys.solver_params.stimLength)
            case false
                disp('No stimulus artifact to be removed ...')
        end
        
    case 'regularization'
        switch sys.solver_params.regularization_flag
            case true
            switch sys.solver_params.lambda
                case 0
                    disp('Maximum Delta F / F reached: started penalization ...')
                otherwise
                    switch sys.solver_params.pen_norm
                        case 'vec_norm'
                            fprintf('Regularizing the elementwise %d norm ... \n',sys.solver_params.dff_qnorm)
                        case 'p_norm'
                            fprintf('Regularizing the %d-%d norm ... \n',sys.solver_params.dff_pnorm,sys.solver_params.dff_qnorm)
                        case 'nuclear'
                            disp('Regularizing the nuclear norm of segment intensities ... \n')
                        otherwise
                            error('penalty type not specified')
                    end
            end
                case false
            disp('No regularization on Delta F / F ...')
        end
    case 'baseline'
        disp('Baseline estimated ...')
        
    case 'ref_image'
        disp('Reference image is not updated ...')
        
    case 'dff break'
        fprintf('Delta F / F changed less than %d percent: Algorithm stops at iteration %d \n',sys.opts.convergence_threshold,sys.control_params.current_iter);
    
    case 'log-likelihood break'
        fprintf('Log likelihood of the estimates changed less than %d percent in the last %d iterations: Algorithm stops at iteration %d \n',sys.opts.convergence_threshold, sys.window_length ,sys.control_params.current_iter);
    
    case 'current frame'
    if mod(sys.control_params.frame,sys.solver_params.verbose) == 0 || sys.control_params.frame == 1 || sys.control_params.frame == sys.opts.T
         fprintf('Writing frame number %d \n',sys.control_params.frame);
    end
    
    case 'objective'
        if sys.control_params.current_iter == 1
            disp('Objective function is calculated without penalty')
            if sys.solver_params.regularization_flag
                disp('Warning: Penalization is also in effect!!')
            end
        end
    
    case 'max_iter'
        if sys.control_params.current_iter == sys.opts.NumIter_dynamic
        disp('Maximum number of iterations reached !!')
        end
    case 'baseline_model_changed'
        if sys.solver_params.regularization_flag
        switch sys.opts.dff_model
            case 'regularized'
            otherwise
            disp('Regularization option for dF/ F is active: Baseline model changed to multiplicative!')
        end
        end
        
    case 'baseline_model'
        error('Invalid baseline model!')
    case 'perSeedBaseline'
        disp(['Per segment baseline estimated at iteration ' num2str(sys.control_params.current_iter)])        
    case 'ReconAlg'
        error('Reconstruction algorithm not defined!')
    case 'NullBaseline'
        disp('Projected data onto the subspace spanned by columns of PS')
    otherwise
    error('Error message type not defined')   
end



end