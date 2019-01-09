function calib = SLAPmi_AI2AO(hSLAPMi)

%ensure timing is updated;
SLAPmi_timing(hSLAPMi);

%Parameters
AOrate = hSLAPMi.galvos.AOrate;
AIrate = AOrate;                 %AIrate must divide AOrate!
NSampsOut = 1e5;
NSampsIn = round(NSampsOut * AIrate/AOrate);
nGalvos = 6; 

import dabs.ni.daqmx.*

galvosOut = most.util.safeCreateTask('galvosOut');
for Gn = 1:nGalvos
    galvosOut.createAOVoltageChan(hSLAPMi.scannerDaqName, 4*(Gn-1), ['galvo' int2str(Gn)], -10, 10);
end 
galvosOut.cfgSampClkTiming(AOrate,'DAQmx_Val_ContSamps');

StartDO = most.util.safeCreateTask('StartDO');
StartDO.createDOChan('DEV6', 'port0/line0' , 'StartDO_TTL');
StartDO.writeDigitalData(false, 1, true);

Vo = -9.6:0.2:9.6;
Vi = nan(6, length(Vo));
galvosOut.cfgOutputBuffer(NSampsOut);
galvosOut.writeAnalogData(Vo(1)*ones(NSampsOut, nGalvos));
start(galvosOut);
disp('Settling...')
pause(3);
stop(galvosOut);
    
for Vix = 1:length(Vo)
    disp(['Voltage:' num2str(Vo(Vix))])
    galvosOut.cfgOutputBuffer(NSampsOut);
    galvosOut.writeAnalogData(Vo(Vix)*ones(NSampsOut, nGalvos));
    f_actual = hSLAPMi.acqGalvoPosSampsAsync(NSampsIn,AIrate,'/DEV5/PFI1');
    start(galvosOut);
    StartDO.writeDigitalData(true, 1, true); %trigger
    pause(NSampsOut/AOrate +0.2);
    stop(galvosOut);
    StartDO.writeDigitalData(false, 1, true); %reset trigger
    
    %read data
    D_raw = hSLAPMi.endGalvoAsyncAcq();
    Vi(:,Vix) = mean(D_raw(ceil(end/2):end,:),1);
end


%Don't make the galvos hang out at +10V
galvosOut.cfgOutputBuffer(NSampsOut);
galvosOut.writeAnalogData(zeros(NSampsOut, nGalvos));
pause(NSampsOut/AOrate +0.2);
stop(galvosOut);

try
    most.idioms.safeDeleteObj([galvosOut, StartDO]);
    hSLAPMi.status = 'Idle';
catch
    keyboard
end

for Gn = 1:nGalvos
    monotone = diff(Vi(Gn,:))>1e-5;
    ix1 = find(~monotone(1:floor(end/2)), 1, 'last')+1;
    if isempty(ix1)
        ix1 = 1;
    end
    monotone(1:floor(end/2)) = true;
    ix2 = find(~monotone, 1, 'first');
    if isempty(ix2)
        ix2 = size(Vi,2);
    end
    select = ix1:ix2;
    hSLAPMi.calib.galvos.AO2AI{Gn} = griddedInterpolant(Vo(select), Vi(Gn,select));
    hSLAPMi.calib.galvos.AI2AO{Gn} = griddedInterpolant(Vi(Gn,select), Vo(select));
end

figure, plot(Vi')

calib = hSLAPMi.calib;