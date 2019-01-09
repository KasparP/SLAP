function [hF, W,H] = SLAPMi_PCA (scandata, k)
    
    y = reshape([scandata.frames.pmtData], size(scandata.frames(1).pmtData,1), size(scandata.frames(1).pmtData,2),[]);
    y = squeeze(y(:,1,:)); %TEMPORARY ONLY ONE CHANNEL
    
    y(isnan(y)) = 0;
    hF = [];
    if nargin<2
        k =3; %number of components to plot
    end

    %perform PCA of scandata and plot
%     [coeff,score,latent] = pca(y);
%     hF = [hF, figure('Name', ['PCA of ' scandata.metadata.galvoDataFileName])];
%     subplot(3,1,1)
%     plot(1e3*scandata.metadata.framePeriod*(1:size(y,2)), coeff(:,1:k)); ylabel('Component Loading'); xlabel('Time (ms)')
%     subplot(3,1,2)
%     plot(score(:,1:k)); ylabel('Component Loading'); xlabel('Lxl')
%     subplot(3,1,3); 
%     plot(latent(1:k)); set(gca,'YScale','log'); ylabel('Explained Variance'); xlabel('Component')
    
    %perform vanilla NMF and plot
    rng(1986) %for reproducible results
    [W,H] = nnmf(y,k); 
    hF = [hF, figure('Name', ['Vanilla NMF of ' scandata.metadata.galvoDataFileName], 'pos', [809 49 784 1068])];
    subplot(2,1,1)
    plot(1e3*scandata.metadata.framePeriod*(1:size(y,2)), H(1,:)); ylabel('Component Loading (offset for visibility)'); xlabel('Time (ms)')
    maxlevel = max(H(1,:));
    hold on
    for fac  = 2:k
        plot(1e3*scandata.metadata.framePeriod*(1:size(y,2)), maxlevel - min(H(fac,:)) + H(fac,:)); 
        maxlevel = maxlevel + max(H(fac,:)) - min(H(fac,:));
    end
    set(gca, 'ytick', [])
    subplot(2,1,2)
    plot(W); ylabel('Component Loading'); xlabel('Lxl')
end