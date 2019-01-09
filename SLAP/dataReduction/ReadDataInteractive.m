import slapmi
[fn dr] = uigetfile('*.gdat');
[pathstr,name,ext] = fileparts([dr filesep fn]);
[metaData, pmtData, galvoData, saveData] = SLAPMi.openDataFiles([pathstr name]);

calib = saveData.calib;
Dalign = pmtData.data(:,2);

%for every ~second of data
G = 2; %A time shift parameter: where to take samples from relative to the peak PMT current
dS = 2e7;  % #samples to process for FFT
fs = 1/pmtData.dt;  %sampling frequency

laserframes = zeros(1,ceil(length(pmtData.data)*5.01e6/fs)); %, 'uint64'); %indexes of laser starts
laserframeix = 0;

n2 = 2^nextpow2(dS);
df = fs / n2; %frequency resolution of fft
F = fft(pmtData.data(1:min(dS, end), 2), n2);

%find maximum near 5MHz
fC = round(5e6/df);
[~,maxind] = max(abs(F(fC+1 + (-1000:1000)))); %+1 because the first bin of F is 0
fActual = 5e6 + df*(maxind-1001);
samplesPerCycle = fs/fActual;

B = repmat(0:samplesPerCycle-1, 10001, 1);

pStarts = round(G+(1:samplesPerCycle:10000*samplesPerCycle));
pStarts = pStarts(pStarts<length(Dalign)-samplesPerCycle);
while ~isempty(pStarts)
        A = repmat(pStarts', 1, size(B,2));
        C = Dalign(A+B(1:length(pStarts),:)); %aligned traces to get phase
        Cavg = median(C,1); %Cavg = trimmean(C,60,1);
        [~,P] = max(diff(Cavg));     
        laserframes(laserframeix + (1:length(pStarts))) = pStarts+P-G;
        laserframeix = laserframeix+length(pStarts);
        
        pStarts = round(double(laserframes(laserframeix)) + samplesPerCycle/2 + (0:samplesPerCycle:10000*samplesPerCycle));
        pStarts = pStarts(pStarts<length(Dalign)-samplesPerCycle);
end

laserframes = laserframes(1:laserframeix);

%%deconvolve
%tau = 1.9;
%D = NN_KP(D, tau);


%Generate Pixel Data
nSamples2sum = 6;
[A, B] = meshgrid(laserframes, 0:nSamples2sum-1);
pixeldata = sum(Dalign(A+B),1);

%generate Galvo positions at the pixel data locations
dG = galvoData.t(end)/length(galvoData.t);
galvoframes = 1+((laserframes-1)*pmtData.dt -galvoData.t(1))/dG;
galvodata =nan(length(pixeldata), 6);
for Gn = 1:6
    tmp = smooth(galvoData.data(:,Gn), 15);
    galvodata(:,Gn) = calib.galvos.AI2AO{Gn}(interp1(tmp, galvoframes));
end

%subdivide data into frames
framePeriod = saveData.framePeriod;
sampsPerFrame = framePeriod/pmtData.dt;
nFrames = floor(length(Dalign)/sampsPerFrame);
pxPerFrame = round(saveData.framePeriod*fActual);
frameStart = 0;
for frame = 1:nFrames
    frameEnd = find(laserframes(frameStart+pxPerFrame+ [-2:2]) <= frame*sampsPerFrame, 1, 'last') + frameStart + pxPerFrame- 3;
    
    frames(nFrames - frame +1).pmtData = pixeldata(frameStart+1:frameEnd); 
    frames(nFrames - frame +1).galvoData = galvodata(frameStart+1:frameEnd,:);
    
    line = saveData.lineIDs(ceil((laserframes(frameStart+1:frameEnd)-(frame-1)*sampsPerFrame)*length(saveData.lineIDs)/sampsPerFrame));
    
    %populate frame data
    frames(nFrames - frame +1).Vx = nan(size(line));
    frames(nFrames - frame +1).Vx(line<=2) = frames(nFrames - frame +1).galvoData(frames(nFrames - frame +1).Vx(line<=2), 1);
    frames(nFrames - frame +1).Vy(line<=2) = frames(nFrames - frame +1).galvoData(frames(nFrames - frame +1).Vy(line<=2), 2);
    frames(nFrames - frame +1).Vx(line>=3) = frames(nFrames - frame +1).galvoData(frames(nFrames - frame +1).Vx(line>=3), 3);
    frames(nFrames - frame +1).Vy(line>=3) = frames(nFrames - frame +1).galvoData(frames(nFrames - frame +1).Vy(line>=3), 4);
    frames(nFrames - frame +1).Z = zeros(size(line));
    frames(nFrames - frame +1).line = line;
    
    frameStart = frameEnd;
end
frames = flipud(frames); %simpler than preallocating a struct






% D = smooth(D,3);
% M = median(D);
% cutoff = M + 0.011; %15 mW above median
% 
% Dcut = D;
% Dcut(Dcut<cutoff) = -10;
% 
% %identify events
% events = findpeaks(Dcut-M);
% 
% 
% p_cut = M+0.025;
% Dcut(Dcut<p_cut) = -10;
% p_events = findpeaks(Dcut-M);
% 
% super_cut = M + 0.240; %superpulse cutoff
% Dcut(Dcut<super_cut) = -10;
% s_events = findpeaks(Dcut-M);
% 
% disp(name)
% disp('# microevents: ')
% disp(length(events))
% 
% disp('# photons: ')
% disp(length(p_events))
% 
% disp('# Superpulses: ')
% disp(length(s_events))
