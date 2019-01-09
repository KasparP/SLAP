function sys = initialize_spikes( sys )
initialization = 'all_ones';
%initialization = 'back_projection';
filtered_data  = imgaussfilt(sys.input.y,[eps 50]);
m = min(filtered_data,[],2);
y = sys.input.y - repmat(m,1,sys.opts.T);

switch initialization
    case 'all_ones'
        X = ones(sys.opts.p,sys.opts.T);
        if isfield(sys.input,'fusedInto')
            ToRemove = sys.input.fusedInto == 0;
            X(ToRemove,:) = 0;
        end
        
        sys.output.spikes = filter(sys.opts.theta,1,X,sys.opts.filter_delays,2);
%         sys.output.spikes = 0.05*ones(sys.opts.p,sys.opts.T);
    case 'back_projection'
        X = ones(sys.opts.p,sys.opts.T);
        spikes1 = filter(sys.opts.theta,1,X,sys.opts.filter_delays,2);
        
        spikes2 = zeros(size(spikes1));
        X = sys.output.PS'*y;
        for p=1:sys.opts.p
            spikes2(p,:) = NND(X(p,:)',10*sys.opts.tau*sys.opts.tau_spacing(1));
        end
        spikes2(:,2:end) = imgaussfilt(spikes2(:,2:end),[eps 2]);
        
        spikes2 = (spikes2-repmat(mean(spikes2,2),1,size(spikes2,2)))./std(spikes2(:,2:end),0,2); %zscore
        spikes2 = max(spikes2-2,0)+2;
        spikes2(:,[1:5 end-10:end]) = 0;
        
        sys.output.spikes = spikes1 + spikes1.*spikes2/2;
        
        %     sys.output.spikes = sys.output.spikes + 0.01*rand(size(sys.output.spikes));
        %         sigma = 2;
        %         x = 1:11;
        %         h = exp(-(x-6).^2/sigma^2/2);
        %         A = filter(h,1,[sys.output.spikes(:,2:end), zeros(size(sys.output.spikes,1),5)] ,[],2);
        %         sys.output.spikes(:,2:end) = A(:,6:end);
    case 'nnls'
        load baseline_nls
        sys.output.spikes = filter(sys.opts.theta,1,X,sys.opts.filter_delays,2);
        sys.output.spikes(sys.output.spikes<0.01) = 0.01;
end
sys.output.F = filter(1,sys.opts.theta,sys.output.spikes,sys.opts.filter_delays,2);

opts.initialization = initialization;
SLAPMi_messages('initialization',opts);

end

