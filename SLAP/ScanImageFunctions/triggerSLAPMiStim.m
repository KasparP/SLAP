function triggerSLAPMiStim(hSLAPMi)
    hTaskTriggerOutOnDmd = most.util.safeCreateTask('TriggerOutOnDemand');
    hTaskTriggerOutOnDmd.createDOChan(hSLAPMi.beamDaqName, 'PFI0' , 'TriggerOutOnDmd');
    hTaskTriggerOutOnDmd.writeDigitalData(false, 1, true);
    hTaskTriggerOutOnDmd.writeDigitalData(true, 1, true);
    disp(['Stimulus triggered at ' datestr(now)])
    hTaskTriggerOutOnDmd.writeDigitalData(false, 1, true); %turn off beam
    delete(hTaskTriggerOutOnDmd);
end