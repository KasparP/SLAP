function calib = SLAPmi_AI2AO_piezo(hSLAPMi)

calib = hSLAPMi.calib;

nV = 20;
Vs = linspace(0,10, nV);

samplerate = 5e5;
nSamps = 5e4;

settletime = 0.5;

hTaskAO = most.util.safeCreateTask('SLAPMIPiezo');
hTaskAO.createAOVoltageChan(hSLAPMi.piezoDaqName, 0);

hTaskAI = most.util.safeCreateTask('SLAPMIPiezoFeedback');
hTaskAI.createAIVoltageChan(hSLAPMi.piezoDaqName, 0);
            
hTaskAO.cfgSampClkTiming(samplerate, 'DAQmx_Val_FiniteSamps', 10);
hTaskAO.cfgOutputBuffer(10); %N_AO_piezo-1 in order to ensure proper retriggering; not tested whether this is necessary

hTaskAI.cfgSampClkTiming(samplerate, 'DAQmx_Val_FiniteSamps', nSamps);
hTaskAI.cfgInputBuffer(nSamps);

Z = nan(length(Vs),1);
for Vix = 1:length(Vs)
    disp(['Position ' int2str(Vix)]);
    hTaskAO.writeAnalogData(Vs(Vix)*ones(10,1));
    hTaskAO.start();
    pause(settletime);
    hTaskAO.stop();
    
    hTaskAI.start();
    pause(nSamps/samplerate);
    Zdata = hTaskAI.readAnalogData(nSamps);
    Z(Vix) = mean(Zdata);
    hTaskAI.stop();
end

%delete Tasks
hTaskAO.abort(); hTaskAI.abort();
most.idioms.safeDeleteObj(hTaskAO);
most.idioms.safeDeleteObj(hTaskAI);

calib.piezo.AO2AI = griddedInterpolant(Vs, Z, 'linear', 'none');
calib.piezo.AI2AO = griddedInterpolant(Z, Vs, 'linear', 'nearest');

figure('name', 'Piezo AO vs AI calibration'), plot(Vs, Z); xlabel('Out voltage'), ylabel ('In voltage')

[FileName,PathName] = uiputfile([hSLAPMi.dataDir '\Calibration\*.cal'], 'Save your calibration file');
if FileName
    save([PathName filesep FileName], 'calib');
end




