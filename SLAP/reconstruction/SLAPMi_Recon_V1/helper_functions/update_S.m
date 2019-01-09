function sys = update_S(sys)
        if sys.control_params.iter_S < sys.solver_params.num_update_S
            if sys.iter_S == 1 && length(indices)~= 1 %size(X,1)
               sys.input.S(sys.input.S_zero_indices) = 1;
            end
            
            sys.input.S(:,indices) = update_S_RL(sys.input.P,sys.input.S,sys.input.y,sys.dff,sys.support_S,indices,sys.NumIter_S);        
            sys.output.PS(:,indices) = sys.input.P*sys.input.S(:,indices);
        end

end