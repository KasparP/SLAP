function [Pup, Pdn, uplines, downlines, Zdes] = linePSF_3D(scandata, optsin)
%Generates a PSF for a drifting Z scan, taking into account that the plane of focus moves over the course of each frame

opts.PSFtype = {'''full''', '''delta'''};  tooltips.PSFtype = 'The shape of the PSF on the diffraction-limited axes';
if nargin>1 %UPDATE WITH USER-SPECIFIED OPTIONS
    for field = fieldnames(optsin)'
         opts.(field{1}) = optsin.(field{1});
    end
else
   opts = optionsGUI(opts, tooltips);
end

framesPerVol = size(scandata.frames,1);
nVol = size(scandata.frames,2);

%find the average Z position of each measurement on the UP and DOWN strokes
%of the piezo
if scandata.metadata.ZBiDi
    uplines = 1:floor(framesPerVol/2);
    downlines = floor(framesPerVol/2)+1:framesPerVol;
else %unidirectional Zscan
    duty = 0.85;
    if isfield(scandata.metadata.piezo, 'duty')
        duty = scandata.metadata.piezo.duty;
    end
    uplines = [(-1:0)+framesPerVol 1:floor(duty*framesPerVol)];
    downlines = [];
end


Pup = getPSF(uplines);
Pdn = [];
if ~isempty(downlines)
    Pdn = getPSF(downlines);
end


function P = getPSF(sel_lines)

Zavg = mean([scandata.frames(sel_lines,:).Z],2)*40; %*40 converts volts to microns
Zavg = Zavg - mean(Zavg);

if strcmpi(opts.PSFtype, 'full')
    psfOpts.Zlxl = Zavg;
    P = linePSF_full(scandata, psfOpts);
elseif strcmpi(opts.PSFtype, 'delta')
    P = linePSF_delta(scandata);
else
    error('invalid PSF type')
end

for line = 1:4
    V = reshape(full(sum(P.P(scandata.line==line,:),1)), length(P.coords{1}), length(P.coords{2}), []);
    figure('Name', ['Line ' int2str(line)]), imshow3D(155*V./max(V(:)))
end

end


end