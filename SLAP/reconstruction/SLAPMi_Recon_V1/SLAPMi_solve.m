function sys = SLAPMi_solve(optsin, fns)
% SLAPMi_solve
%inputs:
%           OPTSIN     settings for the options.
%           FNS        complete paths to problemdata files to solve.
%           if either input is not supplied, a GUI will be generated
%
%output:
%           SYS_RECON  A structure containing the solver output

%User options
opts.tau = 10;              tooltips.tau = 'Decay time constant';
opts.dampar = 1.5;          tooltips.dampar = 'Residuals less than dampar S.D.s will have less effect on updates; typically 0-2';
opts.numIters = 100;        tooltips.numIters = 'Maximum number of iterations';
opts.minNumIters = 50;      tooltips.minNumIters = 'Minimum number of iterations';
opts.darkNoise = 0.01;      tooltips.darkNoise = 'Dark noise, in photons, added to all expected rates';
opts.stimArtifact = 0;      tooltips.stimArtifact = 'to compensate for a stimulus artifact, set to approximate artifact length';
opts.clipGradients = true;  tooltips.clipGradients = 'Clip multiplicative updates to 0.1<X<10 per iteration';
opts.offFocusSpikes = false;tooltips.offFocusSpikes = 'Allow fitting of spikes for scattered excitation light, useful for deep in vivo imaging';
opts.baselineMult = 1;      tooltips.baseline = 'Multiply the estimated baseline by this amount, typically 1 or less';
opts.ReconAlg = 'RL';       tooltips.baseline = 'Reconstruction algorithm: Richardson-Lucy (RL) or Nonnegative Least Squares (nnls)';
opts.correctBleaching = false;        tooltips.correctBleaching = 'fit a rank-1 temporally varying baseline to correct for global bleaching';
opts.threshold_spikes = false; tooltips.threshold_spikes = 'threshold spikes smaller than 1e-10';
opts.clipDFF = inf;          tooltips.clipDFF = 'Clip DFF values above this; set to <=0 or Inf for no maximum';

if nargin  && ~isempty(optsin) %UPDATE WITH USER-SPECIFIED OPTIONS
    for field = fieldnames(optsin)'
        opts.(field{1}) = optsin.(field{1});
    end
else
    opts = optionsGUI(opts, tooltips);
end

%USER SELECT FILE
if nargin<2 || isempty(fns)
    basedir = 'E:\SLAPMiData';
    [filename,directory] = uigetfile([basedir filesep '*PROBLEMDATA*.mat'], 'Select your problem data', 'multiselect', 'on');
    if isempty(filename); return; end
    if ~iscell(filename)
        fns = {[directory filesep filename]};
    else
        fns = strcat(directory, filename);
    end
end
fns = sort_nat(fns);
numFiles = length(fns);
hF = [];
for ix = 1:numFiles
    [directory,filename] = fileparts(fns{ix});
    disp(['Processing: ' filename])
    
    %load problemdata
    load([directory filesep filename])
    if ~(exist('P', 'var')==1)
        disp('Loading saved P matrix...');
        load([directory filesep fnP]);
    end
    
    %INITIALIZE SOLVER
    sys = struct;
    sys.opts = opts;
    sys.opts.filename = [filename '.mat'];
    sys.opts.dr = directory;
    sys.input.P = P.P;
    sys.opts.dim(1) = length(P.coords{1});
    sys.opts.dim(2) = length(P.coords{2});
    sys.opts.Psize = length(P.coords{3});
    %%%%temp
    sys.opts.SLM_minimum = 5e-3;
    S.mask(abs(S.mask - min(S.mask(:)))<0.001) = sys.opts.SLM_minimum;
    if isfield(S,'mask')
        sys.output.SLM_mask = max(S.mask,sys.opts.SLM_minimum);
        temp_SLM = repmat(sys.output.SLM_mask(:),sys.opts.Psize,1);
        A = spdiags(temp_SLM,0,size(sys.input.P,2),size(sys.input.P,2));
        sys.input.P =  sys.input.P * A; %multiply the mask into the projection matrix
    end
    % Define the reference image
    sys.input.ref_image = S.IM;
    sys.input.ref_image(isnan(sys.input.ref_image))=0;
    

    % Append the reference image to segmentation
    if opts.offFocusSpikes
        sys.input.S = [S.seg, ones(size(S.seg,1),1), sys.input.ref_image(:)-sum(S.seg,2)];
    else
        sys.input.S = [S.seg, zeros(size(S.seg,1),1), sys.input.ref_image(:)-sum(S.seg,2)];
%         sys.input.S = [S.seg, sys.input.ref_image(:)];
    end
    if isfield(S,'fusedInto')
        sys.input.fusedInto = S.fusedInto;
    end
        % For Red-Green experiment
    if size(S.data,4) > 1
        A = squeeze(sum(S.data,3));
        sys.input.GreenIm = squeeze(A(:,:,1));
        sys.input.RedIm = squeeze(A(:,:,2));
    end

    %initialize measurements
    y = [scandata.frames.pmtData];
    N = sum(y,2);
    I = isnan(N);
    sys.input.P = sys.input.P(~I,:);
    sys.input.y  = y(~I,:);
    if isfield(opts,'StimOnset')
        sys.opts.StimOnset = opts.StimOnset;
    end
    if opts.stimArtifact %Stimulus Artifact for 2 spot uncaging experiment
        sys.opts.removeStimArtifactFlag =true;
        sys.solver_params.stimLength = opts.stimArtifact;
        [sys.output.stim_artifact,sys.opts.StimOnset] = estimateStimArtifact(scandata, sys.solver_params.stimLength);
        sys.output.stim_artifact = sys.output.stim_artifact(~I,:);
    else
        sys.output.stim_artifact = zeros(size(sys.input.y));
        sys.opts.removeStimArtifactFlag = false;
    end
    sys.output.PS = sys.input.P*sys.input.S;
    
    if any(sys.input.y(:)<0)
        SLAPMi_messages('negative_meas',sys);
        sys.input.y(sys.input.y<0) = 0;
    end
    sys.input.ref = rmfield(S,{'seg','IM','data','mask'});
    sys.input.scandata = rmfield(scandata,{'frames'});
    if numFiles == 1
         clear P S scandata   % clear big variables
    end
    % Apply options
    sys = initialize_solver(sys,opts);
    
    %ensure solver files are added to path
    [dr, ~] = fileparts(which('SLAPMi_solver'));
    addpath(genpath(dr));
    
    %solve
    for iter_dynamic = 1:sys.opts.NumIter_dynamic
    sys.control_params.current_iter = iter_dynamic;
    sys = calc_gradient(sys);
    sys = update_spikes(sys);
    if sys.solver_params.regularization_flag
        sys = update_regularization(sys);
    end
    converged = check_convergence(sys);
    if converged; break ;end
    
    SLAPMi_messages('current iteration',sys);
    SLAPMi_messages('max_iter',sys);
    end
    if isfield(sys.output,'SLM_mask')
        sys.opts.SLM_threshold = prctile(sys.output.SLM_mask(:),5);
        InsideSLM = sys.output.SLM_mask > sys.opts.SLM_threshold;
        nPlanes = size(sys.input.ref_image,3);
        SLM_extended = repmat(InsideSLM,1,1,nPlanes);
        S_temp = sys.input.S(:,[1:end-2,end]).*SLM_extended(:);
        sys.output.OutsideSLMSeeds = ~any(S_temp);
    end
    valid = sys.input.ref.fusedInto>0;
    PS = full(sum(sys.output.PS)');
    sys.output.dPhotons = zeros(size(sys.output.F));
    sys.output.dPhotons(valid,:) = PS(valid).*sys.output.F(sys.input.ref.fusedInto(valid), :);
    % sys = calc_corr(sys);

    sys = rmfield(sys,{'Lp','Lp2','Ln','Ln2','Gp','penalty'});
    recon_conditions(sys); %save the reconstruction
    
    if (isfield(opts,'plot') && opts.plot) % || (~isfield(opts,'plot') && ix<=8)
        try
            close(hF);
        end
        opts.appendName = AppendName(sys);
        opts.MovieType = 'min(dFF,Z)';
        opts.t0 = 508;
        opts.t_end =  1600;
        hF = SLAPMi_Plot(sys,opts);
    end

end
