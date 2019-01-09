function calibrateTimingUsingSLM(hSLAPMi)
%fine calibrate the slapmi timing.
%if this is done correctly, only rigid shifts should be required thereafter
%to align for motion

%There must be a fluorescent sample under the center of the FOV

%write a dot to the center of the SLM
[aa,bb] = meshgrid(1:512);
r = 4;
pattern = sqrt((aa-256.5).^2 + (bb-256.5).^2)<r;

originalRes = hSLAPMi.res;
originalFrames = hSLAPMi.framesToCollect;
originalUser = hSLAPMi.user;

hSLAPMi.user = 'Tests';
hSLAPMi.logData = true;

%take a picture with a very slow scan
resSlow = 9000;
hSLAPMi.res = resSlow;
hSLAPMi.framesToCollect = 100;
hSLAPMi.fileNameStem = ['SLOW' datestr(now, 'DDHHMMSS')];
hF = SLAPmi_cal_Galvos3(hSLAPMi);
waitfor(hF);

SLM_drawPattern(pattern);
hSLAPMi.fileNameApp = 0;
fn1a =  hSLAPMi.galvoDataFileName;
hSLAPMi.initAcq(); pause(0.1);
hSLAPMi.armAcq(); pause(0.5);
hSLAPMi.triggerAcq(); pause(1); %allow acquisition to complete
hSLAPMi.endAcq();
disp('Waiting for writes to complete...'); pause(3);

SLM_drawPattern(false(size(pattern)));
hSLAPMi.fileNameApp = 1;
fn1b =  hSLAPMi.galvoDataFileName;
hSLAPMi.initAcq(); pause(0.1); 
if strcmpi(fn1a, fn1b)
    keyboard
end
hSLAPMi.armAcq(); pause(0.5);
hSLAPMi.triggerAcq(); pause(1); %allow acquisition to complete
hSLAPMi.endAcq();


%take a picture with a very fast scan
resFast = 1300;
hSLAPMi.res = resFast;
hSLAPMi.framesToCollect = 100;
hSLAPMi.fileNameStem = ['FAST' datestr(now, 'DDHHMMSS')];
hF = SLAPmi_cal_Galvos3(hSLAPMi);
waitfor(hF);

SLM_drawPattern(pattern);
hSLAPMi.fileNameApp = 0;
fn2a =  hSLAPMi.galvoDataFileName;
hSLAPMi.initAcq(); pause(0.01);
hSLAPMi.armAcq(); pause(0.5);
hSLAPMi.triggerAcq(); pause(1); %allow acquisition to complete
hSLAPMi.endAcq();
disp('Waiting for writes to complete...'); pause(2)
SLM_drawPattern(false(size(pattern)));
hSLAPMi.fileNameApp = 1;
fn2b =  hSLAPMi.galvoDataFileName;
if strcmpi(fn2a, fn2b)
    keyboard
end
hSLAPMi.initAcq(); pause(0.01);
hSLAPMi.armAcq(); pause(0.5);
hSLAPMi.triggerAcq(); pause(1); %allow acquisition to complete
hSLAPMi.endAcq();

%reduce scans
disp('Waiting for writes to complete...'); pause(2)
opts.path = fn1a;
SD1a = SLAPMi_reduce(opts);
opts.path = fn1b;
SD1b = SLAPMi_reduce(opts);
opts.path = fn2a;
SD2a = SLAPMi_reduce(opts);
opts.path = fn2b;
SD2b = SLAPMi_reduce(opts);

%find the center of the aperture under each condition
y1 = nanmean([SD1a.frames.pmtData],2) - nanmean([SD1b.frames.pmtData],2);
y2 = nanmean([SD2a.frames.pmtData],2) - nanmean([SD2b.frames.pmtData],2);
delay = nan(1,4);
for line = 1:4
    L1 = y1(SD1a.line==line); L1([1:100, end-99:end]) = 0;
    V1 = SD1a.Vx(SD1a.line==line);
    [~, center1] = max(smooth(L1, floor(resSlow/20)));
    L2 = y2(SD2a.line==line); L2([1:100, end-99:end]) = 0;
    V2 = SD2a.Vx(SD2a.line==line);
    [~, center2] = max(smooth(L2, floor(resFast/20)));
    
    Gcenter = V1(center1);
    center2_exp = interp1(V2, 1:length(V2), Gcenter);
    delay(line) = (center2-center2_exp)/SD1a.laserFreq;
end

%update the galvo AIDelay
newdelay  = hSLAPMi.calib.galvos.AIdelay + [mean(delay(1:2)) mean(delay(1:2)) mean(delay(3:4)) mean(delay(3:4))];
keyboard

hSLAPMi.calib.galvos.AIdelay = newdelay;
hSLAPMi.res = originalRes;
hSLAPMi.framesToCollect = originalFrames;
hSLAPMi.user = originalUser;
