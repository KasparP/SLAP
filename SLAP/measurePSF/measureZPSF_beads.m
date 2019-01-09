function measureZPSF_beads (hSICtl, hSLAPMi) 
%take a series of acquisitions while stepping through Z
%used to obtain data for calculating the 3D PSF

Zscandir = 'E:\SLAPMiData\PSF\Zscans';
if ~exist(Zscandir, 'dir')
    mkdir(Zscandir);
end
cd(Zscandir);

Zcenter = hSICtl.hModel.hFastZ.positionTarget; %hSICtl.hModel.hMotors.motorPosition(4);
Zlist = Zcenter + (-6:0.2:6);
if any(Zlist<0)
    error('Move the piezo to a center position to obtain a centered stack')
end

hSLAPMi.framesToCollect = 100;

hSLAPMi.logData = true;
hSLAPMi.initAcq()
hSLAPMi.armAcq();
try
    for Z_ix = 1:length(Zlist)
        hSICtl.hModel.hFastZ.positionTarget = Zlist(Z_ix);
        hSICtl.changedMotorPosition;
        pause(0.3) %allow motor to settle
        
        hSLAPMi.triggerAcq();
        pause(0.3) %allow acquisition to complete
        %hSLAPMi.abortAcq(obj)
    end
    hSLAPMi.endAcq();
    
catch ME
    hSLAPMi.abortAcq();
    hSLAPMi.endAcq();
    keyboard
end

