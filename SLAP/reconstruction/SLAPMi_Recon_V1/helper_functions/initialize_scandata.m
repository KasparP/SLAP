sys.opts.filename = filename;
sys.opts.dr = directory;

%%
middle = ceil(size(S.IM,3)/2);
sys.opts.dim(1) = length(P.coords{1});
sys.opts.dim(2) = length(P.coords{2});
sys.opts.Psize = length(P.coords{3});

%% Projection Matrix
sys.input.P = P.P;
sys.opts.SLM_minimum = 5e-3;
if isfield(S,'mask')
    sys.output.SLM_mask = max(S.mask,sys.opts.SLM_minimum);
    temp_SLM = repmat(sys.output.SLM_mask(:),sys.opts.Psize,1);
    A = spdiags(temp_SLM,0,size(sys.P,2),size(sys.P,2));
    sys.input.P =  sys.input.P * A;
end

%% Define the reference image
sys.input.ref_image = S.IM;
sys.input.ref_image(isnan(sys.input.ref_image))=0;


%% Append the reference image
sys.S = [S.seg, ones(size(S.seg,1),1), sys.input.ref_image(:)-sum(S.seg,2)];

%% Measurements

y = [scandata.frames.pmtData];

T = size(y,2);
N = sum(y,2);
I = isnan(N);
sys.input.P = sys.input.P(~I,:);
sys.input.y  = y(~I,:);

if sys.opts.removeStimArtifactFlag
    sys.output.stim_artifact = estimateStimArtifact(scandata, sys.solver_params.stimLength);
    sys.output.stim_artifact = sys.output.stim_artifact(~I,:);
end

if any(sys.input.y(:)<0)
    SLAPMi_messages('negative_meas',sys);
end
sys.input.y(sys.input.y<0) = 0;

%%
clear P S S0 scandata