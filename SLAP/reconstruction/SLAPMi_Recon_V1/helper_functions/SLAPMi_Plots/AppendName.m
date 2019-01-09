function name = AppendName(sys)
    name = [];
    if isfield(sys.opts,'offFocusSpike')
        if sys.opts.offFocusSpike
            name = [name '_' 'offFocusSpike'];
        end
    end
    name = [name '_' 'dampar' num2str(10*sys.solver_params.dampar)];
    name = [name '_' 'tau' num2str(sys.opts.tau)];
    name = [name '_' 'baseline' num2str(100*sys.opts.baselineMult)];
    name = [name '_' num2str(sys.opts.NumIter_dynamic) 'iters'];
    if isfield(sys.opts,'clipDFF') && ~isinf(sys.opts.clipDFF)
    name = [name '_clipDFF' num2str(sys.opts.clipDFF)];
    end
end