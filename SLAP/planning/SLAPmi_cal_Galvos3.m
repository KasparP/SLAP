function hF = SLAPmi_cal_Galvos3 (hSLAPMi)
if ~isempty(hSLAPMi) && ~strcmpi(hSLAPMi.status, 'idle')
    disp('SLAPmi must be idle to calibrate galvos')
    return
else
    hSLAPMi.status = 'Calibrating Galvos';
end

%ensure timing is updated;
SLAPmi_timing(hSLAPMi);

basedir = hSLAPMi.dataDir;
RUNNING = false;

%Parameters
maxiter = 40;
pausetime = 1; %min(2, 500/freq + 0.1);     %how long to run the galvos, in seconds
freq =  1/hSLAPMi.framePeriod;
nrepeats = max(3, floor((pausetime-0.1)*freq));
AOrate = hSLAPMi.galvos.AOrate;
AIrate = AOrate;                 %AIrate must divide AOrate!
%powerRate = hSLAPMi.sampleRate;
NSampsOut = 4*round(AOrate/(4*freq));
if abs(AOrate/(4*freq) - round(AOrate/(4*freq)))>1e-5
    keyboard
end
NSampsIn = round(NSampsOut * AIrate/AOrate);
nGalvos = 6;
colors = [1 0 0; 0 0 1 ; 0 0.5 0; 1 0 1; 0.5 0.5 0; 0 1 1];

disp(['Actual frequency: ' num2str(AOrate/NSampsOut)]) %display actual frequency

%draw GUI
hF = figure('name', 'Galvo Calibration', 'pos', [960 570 570 540], 'closerequestfcn', @closeFcn, 'dockcontrols', 'off', 'numbertitle', 'off');
hAx1 = axes('parent', hF, 'pos', [0.05 0.05 0.4 0.4]);
hAx4 = axes('parent', hF, 'pos', [0.55 0.05 0.4 0.4]);
hAx5 = axes('parent', hF, 'pos', [0.05 0.55 0.4 0.4]);

hAbort = uicontrol(hF, 'units', 'norm', 'String', 'START');  %restart/abort button
set(hAbort, 'pos', [0.65 0.85 0.2 0.08], 'visible', 'On', 'Callback', @startCB)

hReset = uicontrol(hF, 'units', 'norm', 'String', 'Reset');  %restart/abort button
set(hReset, 'pos', [0.65 0.75 0.2 0.08], 'visible', 'On', 'Callback', @resetCB)

hAccept = uicontrol(hF, 'units', 'norm', 'String', 'Accept');  %restart/abort button
set(hAccept, 'pos', [0.65 0.65 0.2 0.08], 'visible', 'On', 'Callback', @acceptCB)

apCalib = hSLAPMi.calib;
apCalib.galvos.linelength = hSLAPMi.PSF.aperture.linelength * hSLAPMi.aperture/0.4; %The aperture linelength is always calibrated at a size of 0.4
for line = 1:4
    lname = (['line' int2str(line)]);
    apCalib.galvos.offset.(lname).X = hSLAPMi.PSF.aperture.(lname).offsetX;
    apCalib.galvos.offset.(lname).Y = hSLAPMi.PSF.aperture.(lname).offsetY;
end
[desired, hSLAPMi.lineIDs] = galvo_trace_SLAPmi(NSampsOut,hSLAPMi.tiling, apCalib);

plot(hAx1, desired);

x1 = linspace(0,1, NSampsOut+1);
x2 = linspace(0,1, NSampsIn+1);
des_resamp = interp1(x1(1:end-1), desired, x2(1:end-1));
des_x = nan(size(des_resamp));
for Gn = 1:nGalvos
    des_x(:,Gn) = hSLAPMi.calib.galvos.AO2AI{Gn}(des_resamp(:,Gn));
end

%initialize AO; if we've already calibrated before, use that trace
hSLAPMi.galvos.desired = desired;
AO = desired;
% target = des_x - mean(des_x,1);
% offset = mean(des_x,1);
beamAO = zeros(size(AO,1),1);
ID = DataHash({AOrate, desired, 3}); %Hash of the galvo signal
if exist([basedir filesep 'GalvoTraces\' ID '.AO'], 'file')
    load([basedir filesep 'GalvoTraces\' ID '.AO'], '-mat')
end

%SETUP DAQ
import dabs.ni.daqmx.*
galvosOut = most.util.safeCreateTask('galvosOut');
galvosOut.createAOVoltageChan(hSLAPMi.scannerDaqName, hSLAPMi.scannerDaqChannels);
galvosOut.cfgSampClkTiming(AOrate,'DAQmx_Val_ContSamps');
set(galvosOut,'sampClkTimebaseSrc', 'PXI_Trig7');
set(galvosOut,'sampClkTimebaseRate', 5000000);
StartDO = most.util.safeCreateTask('StartDO');
StartDO.createDOChan('DEV6', 'port0/line0' , 'StartDO_TTL');
StartDO.writeDigitalData(false, 1, true);
galvosOut.cfgDigEdgeStartTrig('/DEV5/PFI1');
powerOnDemand = most.util.safeCreateTask('SLAPMIPowerBeamOnDemand');
powerOnDemand.createAOVoltageChan(hSLAPMi.beamDaqName, 0);
powerOut = most.util.safeCreateTask('SLAPMIPowerBeam');
powerOut.createAOVoltageChan(hSLAPMi.beamDaqName, 0);

%lines are at the following phase:
lin_ixs{1} = repmat([true(1, NSampsIn/(4*hSLAPMi.tiling)) false(1, NSampsIn/(4*hSLAPMi.tiling))], 1, 2*hSLAPMi.tiling); %for A and B (and E) galvos
lin_ixs{2} = ~lin_ixs{1};  %for C and D (and F) galvos
des_line = [des_x(lin_ixs{1}, 1:2) des_x(lin_ixs{2}, 3:4) des_x(lin_ixs{1}, 5) des_x(lin_ixs{2}, 6) ];

%extize
errors = nan(nGalvos,maxiter);
iter = -1; %negative number means iterations to run before measuring errors;

abort = false;
p = []; pVec = zeros(1,6);
h = []; %H = [];
aa = [];
bb = [];
KK = cell(1,nGalvos);
lambda = 300;
    function startCB(varargin)
        if RUNNING
            abort = true;
            set(hAbort, 'string', 'START', 'UserData', true)
            RUNNING = false;
        else
            abort = false;
            RUNNING = true;
            set(hAbort, 'string', 'STOP', 'Userdata', false)
            
            %figure out phase delay
            if isempty(p)
                disp('Measuring Impulse Response...')
                getImpulseResponse;
                disp('Allocating convolution matrices...')
                N = size(aa,1);
                DD = convmtx([1 -2 1], N);
                DD(:, 1:2) = DD(:,1:2) + DD(:, end-1:end);
                DD = DD(:,1:N);
                DD = DD'*DD;
                for Gn = 1:nGalvos
                    K = convmtx(h(:,Gn), N); 
                    K = K(1:N,:) + K([N+1:end 1], :);
                    KK{Gn} = pinv(K'*K + lambda*DD)*K';
                end
            end
            
            if isempty(AO)
                AO = circshift(desired, [-p 0]);
            end

            while ~abort && ~get(hAbort, 'UserData')
                iter = iter+1;
                if iter>maxiter
                    abort = true;
                end
                
                AO = min(max(AO,-10),10);
                
                galvosOut.cfgOutputBuffer(NSampsOut);
                galvosOut.writeAnalogData([AO beamAO]);
                f_actual = hSLAPMi.acqGalvoPosSampsAsync(NSampsIn*nrepeats,AIrate,'/DEV5/PFI1');  %#ok<NASGU> %output is necessary!
                start(galvosOut);
                StartDO.writeDigitalData(true, 1, true); %trigger
                pause(pausetime);
                stop(galvosOut);
                StartDO.writeDigitalData(false, 1, true); %reset trigger
                
                %read data
                D_raw = hSLAPMi.endGalvoAsyncAcq();
                
                if abort || get(hAbort, 'UserData')  || iter<1
                    continue
                end
                
                D = zeros(NSampsIn, nGalvos);
                for Gn = 1:nGalvos
                    D(:,Gn) = squeeze(mean(interp1(D_raw(:,Gn), aa+bb),2));
                end
                
                %compute errors
                E{1} = des_line(:,1) -D(lin_ixs{1},1);
                E{2} = des_line(:,2) -D(lin_ixs{1},2);
                E{3} = des_line(:,3) -D(lin_ixs{2},3);
                E{4} = des_line(:,4) -D(lin_ixs{2},4);
                E{5} = des_line(:,5) -D(lin_ixs{1},5);
                E{6} = des_line(:,6) -D(lin_ixs{2},6);
                errors(1:2,iter) = sqrt(mean((des_line(:,1:2) -D(lin_ixs{1},1:2)).^2));
                errors(3:4,iter) = sqrt(mean((des_line(:,3:4) -D(lin_ixs{2},3:4)).^2));
                errors(5,iter) = sqrt(mean((des_line(:,5) -D(lin_ixs{1},5)).^2));
                errors(6,iter) = sqrt(mean((des_line(:,6) -D(lin_ixs{2},6)).^2));
                
                %UPDATE DRIVE SIGNAL
                alpha = 0.9;
                newAO = AO;
                
                for Gn = 1:nGalvos
                    err = D(:,Gn) - des_x(:,Gn); 
                    err_m = mean(err);
                    err = err-err_m;
                    
                    err_i = KK{Gn} * err;
                    err_i = err_i - mean(err_i);
                    newAO(:,Gn) = AO(:,Gn) - alpha*(err_i(1:size(AO,1))) + 2*mean(E{Gn}); %really should be (AI2AO(err_m) - AI2AO(0))
                end
                
                AO = newAO;
                
                %plotting
                cla(hAx1), cla(hAx4), cla(hAx5)
                for G = 1:nGalvos
                    plot(hAx4, AO(:,G), 'color', colors(G,:))
                    hold(hAx4, 'on')
                    plot(hAx5, errors(G,:), 'color', colors(G,:), 'linewidth', 2); 
                    hold(hAx5, 'on'), scatter(hAx5, iter, errors(G,end), 'MarkerEdgeColor', colors(G,:));
                    plot(hAx1,des_x(:,G),  'k'), hold(hAx1, 'on'), plot(hAx1, D(:,G), 'color', colors(G,:), 'linewidth', 2)
                end
                set(hAx5, 'YScale', 'log'); xlabel(hAx1, 'Position'); xlabel(hAx4, 'Command Signal'); xlabel(hAx5, 'Error History');
                drawnow;
            end
            %aborted
            RUNNING = false;
        end
    end

    function acceptCB(varargin)
        if ~RUNNING
            hSLAPMi.galvos.calFramePeriod = hSLAPMi.framePeriod;
            hSLAPMi.galvos.calPos = [hSLAPMi.scanCenterPoint hSLAPMi.tiling hSLAPMi.aperture];
            hSLAPMi.framePeriod = hSLAPMi.framePeriod;
            hSLAPMi.galvoWaveform = AO;
            [beamWaveform, hSLAPMi.powerWaveformRaw] = pockels_trace_SLAPmi(NSampsOut, hSLAPMi);
            hSLAPMi.beamWaveform = circshift(beamWaveform, -[round(AOrate*hSLAPMi.calib.pockels.wireDelayFast) 0]); %shift the AO to compensate for wire transmission delay
            
            save([basedir filesep 'GalvoTraces\' ID '.AO'], 'AO'); %save the calibrated AO under the hash
            closeFcn();
        end
    end

    function closeFcn (varargin)
        try
            most.idioms.safeDeleteObj([galvosOut, StartDO, powerOnDemand, powerOut]);
            hSLAPMi.status = 'Idle';
        catch
        end
        delete(hF)
    end

    function resetCB (varargin)
        if ~RUNNING
            %extize
            errors = nan(nGalvos,maxiter);
            iter = -1;
            AO = desired;
            abort = false;
            p = [];
            
            cla(hAx1); cla(hAx4) ; cla(hAx5)
        end
    end

    function getImpulseResponse
        dt = 100;
        X = zeros(NSampsOut, 7);
        X(dt,1:6) = 1;
        
        repeat = 2;
        while repeat
            galvosOut.cfgOutputBuffer(NSampsOut);
            galvosOut.writeAnalogData(X);
            fs = hSLAPMi.acqGalvoPosSampsAsync(NSampsIn*nrepeats,AIrate,'/DEV5/PFI1'); %#ok<NASGU>
            
            start(galvosOut);
            StartDO.writeDigitalData(true, 1, true); %trigger
            
            pause(pausetime);
            
            StartDO.writeDigitalData(false, 1, true); %reset trigger
            stop(galvosOut);
            
            D_raw = hSLAPMi.endGalvoAsyncAcq();
            
            %The input and output clocks are not synchronized, we need to
            %resample the input to avoid phase precession
            %for some reason periodogram has issues here??
            
            %find the peak of the first and last impulses
            [~, max1] = max(smooth(D_raw(3*NSampsIn + (1:NSampsIn),1),31)); 
            [~, max2] = max(smooth(D_raw(end-NSampsIn+1:end,1),31));
            max2 = max2+size(D_raw,1)-NSampsIn;
            max1 = max1+3*NSampsIn;
            Ncycles = round((max2-max1)/NSampsIn); %nrepeats; 
            cycleframes = (max2-max1)/Ncycles;
            delaySamps = hSLAPMi.calib.galvos.AIdelay*AIrate;
            [aa,bb] = meshgrid(cycleframes+delaySamps:cycleframes:size(D_raw,1)-NSampsIn, 1:NSampsIn);
            
            D1 = zeros(NSampsIn, nGalvos);
            for Gn = 1:nGalvos
                D1(:,Gn) = squeeze(mean(interp1(D_raw(:,Gn), aa+bb),2));
                h1  = D1(:,Gn) - mean(D1(:,Gn));
                h(:,Gn) = circshift(h1, [floor(-dt*AIrate/AOrate), 0]);
                [~, pVec(:,Gn)] = max(h(:,Gn));
            end
            
            p = round(mean(pVec(1:4)));
            
            if repeat==1 && (any(pVec<AIrate*5e-5) || any (pVec>AIrate*25e-5))
                repeat = strcmpi('Yes', questdlg({'Phase Delay came out strange!'; int2str(pVec) ;'Ensure your inputs are connected. Retry?'}, 'Retry?', 'Yes', 'No', 'Yes'));
                if ~repeat
                    keyboard
                end
            else
                repeat = repeat-1;
            end
        end
    end
end

function [T,Td] = V2T(Vx,Vy,Xoff,Yoff,angle)
T = (Vx - Xoff)*cos(angle) - (Vy-Yoff)*sin(angle);
Td = (Vx - Xoff)*sin(angle) + (Vy-Yoff)*cos(angle);
end

function [Vx,Vy] = T2V(T,Td,Xoff,Yoff,angle) %#ok<DEFNU>
Vx = T*cos(angle) + Td*sin(angle) + Xoff;
Vy = -T*sin(angle) + Td*cos(angle) + Yoff;
end



