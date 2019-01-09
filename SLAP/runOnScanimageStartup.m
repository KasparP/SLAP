function runOnScanimageStartup(varargin)
hSICtl = evalin('base', 'hSICtl');
hSI = evalin('base', 'hSI');

    %turn on SLM control
    SLMcontrol; 
    
    %move piezo to center
    
    hSICtl.hModel.hFastZ.positionTarget = 200;
    hSICtl.changedMotorPosition;
    
    %turn off scanimage averaging
    hSI.hScan2D.hAcq.disableMatlabAveraging = false;
    hSI.hScan2D.hAcq.disableFpgaAveraging = false;%hSI.hScan2D.hAcq.disableFpgaAveraging = true;

end