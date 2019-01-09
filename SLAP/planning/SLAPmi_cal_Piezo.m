function SLAPmi_cal_Piezo (hSLAPMi)
if ~isempty(hSLAPMi) && ~strcmpi(hSLAPMi.status, 'idle')
    disp('SLAPmi must be idle to calibrate Piezo')
    return
else
    hSLAPMi.status = 'Calibrating Piezo';
end

%ensure timing is updated;
SLAPmi_timing(hSLAPMi);

RUNNING = false;

%Parameters
maxiter = 400;
nrepeats = 10;
freq =  1/hSLAPMi.framePeriod;
AOrate = hSLAPMi.galvos.AOrate;
AIrate = AOrate;                 %AIrate must divide AOrate!
NSampsOut = 4*round(AOrate/(4*freq));
if AOrate/(4*freq) ~= round(AOrate/(4*freq))
    keyboard
end
NSampsIn = round(NSampsOut * AIrate/AOrate);
S = 12e-6*AOrate;               %how many samples to shift the changes in time. 5-10 is good, depends on samplerate
A = 1.5;                        %how much to amplify changes in the linear region
S2 = round(NSampsIn/5);         %which changes to amplify
pausetime = 400/freq + 0.1;     %how long to run the galvos
nGalvos = 2; 
colors = parula(nGalvos);


disp(['Actual frequency: ' num2str(AOrate/NSampsOut)]) %display actual frequency


%draw GUI
hF = figure('name', 'Galvo Calibration GUI', 'pos', [809 49 784 1068], 'resize', 'off', 'windowstyle', 'modal', 'closerequestfcn', @closeFcn);
hAx1 = axes('parent', hF, 'pos', [0.05 0.05 0.4 0.25]);   
hAx2 = axes('parent', hF, 'pos', [0.05 0.38 0.4 0.25]);   
hAx3 = axes('parent', hF, 'pos', [0.55 0.05 0.4 0.25]);   
hAx4 = axes('parent', hF, 'pos', [0.55 0.38 0.4 0.25]);
hAx5 = axes('parent', hF, 'pos', [0.05 0.7 0.4 0.25]);  

hAbort = uicontrol(hF, 'units', 'norm', 'String', 'START');  %restart/abort button
set(hAbort, 'pos', [0.65 0.8 0.2 0.03], 'visible', 'On', 'Callback', @startCB)

hReset = uicontrol(hF, 'units', 'norm', 'String', 'Reset');  %restart/abort button
set(hReset, 'pos', [0.65 0.75 0.2 0.03], 'visible', 'On', 'Callback', @resetCB)

hAccept = uicontrol(hF, 'units', 'norm', 'String', 'Accept');  %restart/abort button
set(hAccept, 'pos', [0.65 0.7 0.2 0.03], 'visible', 'On', 'Callback', @acceptCB)


desired = galvo_trace_SLAPmi(NSampsOut,hSLAPMi.calib);  %galvo_trace slapmi produces traces such that the total length (x+y) of each linear segment is 1
plot(hAx1, desired(:, 1:nGalvos));

%SETUP DAQ
import dabs.ni.daqmx.*
galvosOut = most.util.safeCreateTask('galvosOut');
for Gn = 1:nGalvos
    galvosOut.createAOVoltageChan(hSLAPMi.scannerDaqName, 4*(Gn-1), ['galvo' int2str(Gn)], -10, 10);
end 
galvosOut.cfgSampClkTiming(AOrate,'DAQmx_Val_ContSamps');
StartDO = most.util.safeCreateTask('StartDO');
StartDO.createDOChan('DEV6', 'port0/line0' , 'StartDO_TTL');
StartDO.writeDigitalData(false, 1, true);
galvosOut.cfgDigEdgeStartTrig('/DEV5/PFI1');

x1 = linspace(0,2,2*NSampsOut+1);
x2 = linspace(0, 1, NSampsIn+1);
des_resamp = interp1(x1(1:end-1), repmat(desired,2,1), x2(1:end-1));
des_x = des_resamp/4; %the MOD and MON voltages seem to have different gains
des_v = diff(des_x([1:end 1],:),1,1);
des_a = diff(des_v([1:end 1],:),1,1);

%lines are at the following phase:
lin_ixs{1} = [1:(NSampsIn/4)  (NSampsIn/2)+(1:(NSampsIn/4))];
lin_ixs{2} = [(NSampsIn/4)+1:(NSampsIn/2) (3*NSampsIn/4)+1:NSampsIn];
des_line = [des_x(lin_ixs{1}, 1:2) des_x(lin_ixs{2}, 3:4)];

%extize
TDold = zeros(nGalvos, NSampsOut);
TUold = zeros(nGalvos, NSampsOut);
errors = nan(4,maxiter);
iter = 0;
AO = desired;
abort = false;
alphas = ones(1,nGalvos);
p = [];


function startCB(varargin)
    if RUNNING
       abort = true;
       set(hAbort, 'string', 'START', 'UserData', true)
       RUNNING = false;
    else
        RUNNING = true;
        set(hAbort, 'string', 'STOP', 'Userdata', false)
        
        %figure out phase delay
        if isempty(p)
            X = linspace(0,2*pi, NSampsOut+1); X = X(2:end);
            X = [sin(X)'  sin(X+pi)'  sin(X+pi/2)'  sin(X+3*pi/2)'];
            %X = X*0.8;
            %X = X.*repmat(max(AO,[],1) - min(AO,[],1), size(X,1),1);
            
            galvosOut.cfgOutputBuffer(NSampsOut);
            galvosOut.writeAnalogData(X);
            f_actual = hSLAPMi.acqGalvoPosSampsAsync(NSampsIn*nrepeats,AIrate,'/DEV5/PFI1');

            start(galvosOut);
            StartDO.writeDigitalData(true, 1, true); %trigger
            
            pause(pausetime);
            %pause(10)
            
            StartDO.writeDigitalData(false, 1, true); %reset trigger
            stop(galvosOut);
            
            D_raw = hSLAPMi.endGalvoAsyncAcq();
            
            D1 = squeeze(mean(reshape(D_raw(NSampsIn+1:end,:), NSampsIn, nrepeats-1,4),2));
            
            %DFT the input signal and the measured signal to recover phase
            Pdata = sum(D1.*repmat(exp(-1i*2*pi*(1:NSampsOut)/NSampsOut),4,1)')*2/NSampsOut;
            PX = sum(X.*repmat(exp(-1i*2*pi*(1:NSampsOut)/NSampsOut),4,1)')*2/NSampsOut;
            P = angle(Pdata) - angle(PX); %phase difference
            
            p = round(mean(mod(P, 2*pi)) * NSampsOut/(2*pi));
            
            %identify the phase of the sinewave for each channel, relative
            %to the input
            p = 84;
        end
        
        while ~abort && ~get(hAbort, 'UserData')
            iter = iter+1;
            AO = min(max(AO,-9.9),9.9);
            if iter>maxiter
                abort = true;
            end
            if any(AO(:)>=9.5) || any (AO(:)<=-9.5)
                warning('Very large signals being sent to galvos!');
            end
            
            galvosOut.cfgOutputBuffer(NSampsOut);
            galvosOut.writeAnalogData(AO);
            
            f_actual = hSLAPMi.acqGalvoPosSampsAsync(NSampsIn*nrepeats,AIrate,'/DEV5/PFI1');
            start(galvosOut);
            StartDO.writeDigitalData(true, 1, true); %trigger

            pause(pausetime);

            stop(galvosOut);
            StartDO.writeDigitalData(false, 1, true); %reset trigger

            %read data
            D_raw = hSLAPMi.endGalvoAsyncAcq();
           
            if abort || get(hAbort, 'UserData') 
                continue
            end
            
            keyboard; 
            
            D = squeeze(mean(reshape(D_raw(NSampsIn+1:end,:), NSampsIn, nrepeats-1,4),2));
            errors(1:2,iter) = sqrt(mean((des_line(:,1:2) -D(mod(lin_ixs{1}+p-1, NSampsIn)+1,1:2)).^2));
            errors(3:4,iter) = sqrt(mean((des_line(:,1:2) -D(mod(lin_ixs{2}+p-1, NSampsIn)+1,3:4)).^2));
            
            D_align = circshift(D, [-p 0]);
            
            v = nan(NSampsIn, 4); %measured velocity of mirrors, trend filtered
            for G = 1:4
                tmp = l1tf(diff(repmat(D_align(:,G),[3,1])),0.05); %trend filtering
                v(:,G) = tmp(NSampsIn+1:2*NSampsIn);
                keyboard
            end

            %UPDATE DRIVE SIGNAL

            for G = 1:2
                galvoset = ceil(G/2.1);
                
                errslope = 0;
                if iter>5
                    errslope = polyfit(1:5, log(errors(G, iter-4:iter)),1);
                end
                alphas(G) = max(0,min(0.05, errors(G,iter)*(1-(30*errslope(1)))));
                
                P = des_x(:,G) - D_align(:,G);
                D = des_v(:,G) -  v(:,G);
                
                targetUP = double(P>0 & D>0); % & a<des_a'; %where the actual trace is below the desired, the slope is too low, and  %min(2,(0.5+5*P)).*
                targetUP(mod(lin_ixs{galvoset}-S2, NSampsIn)+1) = targetUP(mod(lin_ixs{galvoset}-S2, NSampsIn)+1)*A;
                
                targetUP = interp1(linspace(0,1,length(targetUP)+1), targetUP([1:end,1]), linspace(0,1,NSampsOut+1));
                targetUP = circshift(targetUP(1:end-1), [0 S]);  %empirical shift to accomodate for delay in read/write
                
                targetUP = targetUP + 0.4*(targetUP.*TUold(G,:)); %momentum in conserved change regions
                TUold(G,:) = targetUP;
                
                targetDOWN = double(P<0 & D<0); % & a>des_a'; %pos is high, vel is high, accel is high  %min(2,(0.5-5*P)).*
                targetDOWN(mod(lin_ixs{1}-S2, NSampsIn)+1) = targetDOWN(mod(lin_ixs{1}-S2, NSampsIn)+1)*A;
                
                targetDOWN = interp1(linspace(0,1,length(targetDOWN)+1), targetDOWN([1:end,1]), linspace(0,1,NSampsOut+1));
                targetDOWN = circshift(targetDOWN(1:end-1), [0 S]); %empirical shift to accomodate for delay in read/write
                
                targetDOWN = targetDOWN + 0.4*(targetDOWN.*TDold(G,:)); %momentum in conserved change regions
                TDold(G,:) = targetDOWN;
                
                AO(:,G) = AO(:,G)+alphas(G)*targetUP';
                AO(:,G) = AO(:,G)-alphas(G)*targetDOWN';
                
                
                cla(hAx1), cla(hAx2), cla(hAx3), cla(hAx4), cla(hAx5)
                plot(hAx1,des_x,  'b'), hold(hAx1, 'on'), plot(hAx1, D_align, 'r')
                plot(hAx2,des_a), hold(hAx2, 'on'), plot(hAx2, diff(v([1:end,1],:),1,1), 'r')
                plot(hAx3, targetUP+(G-1)*2, 'k')
                hold(hAx3, 'on')
                plot(hAx3, targetDOWN+(G-1)*2, 'r')
                plot(hAx4, AO(:,G))
                hold(hAx4, 'on')
                plot(hAx5, errors(G,:), 'color', colors(G,:)); hold on, scatter(iter, errors(G,end), 'markeredge', colors(G,:));
                hold(hAx5, 'on')
            end
            set(hAx5, 'YScale', 'log')
            xlabel(hAx1, 'Position'); xlabel(hAx2, 'Acceleration'); xlabel(hAx3, 'Change Sites'); xlabel(hAx4, 'Command Signal'); xlabel(hAx5, 'Error History');

        end
        RUNNING = false;
    end

end

    function acceptCB(varargin)
        hSLAPMi.galvos.calFramePeriod = hSLAPMi.framePeriod;
        hSLAPMi.framePeriod = hSLAPMi.framePeriod;
        hSLAPMi.galvoWaveform = AO;
        hSLAPMi.beamWaveform = pockels_trace_SLAPmi(NSampsOut, hSLAPMi.calib);
        closeFcn();
    end

    function closeFcn (varargin)
        most.idioms.safeDeleteObj([galvosOut, StartDO]);
        try
            hSLAPMi.status = 'Idle';
        catch
        end
        delete(hF)
    end

    function resetCB (varargin)
        if ~RUNNING
            %extize
            TDold = zeros(nGalvos, NSampsOut);
            TUold = zeros(nGalvos, NSampsOut);
            errors = nan(4,maxiter);
            iter = 0;
            AO = desired;
            abort = false;
            alphas = ones(1,nGalvos);
            p = [];
        end
    end

end