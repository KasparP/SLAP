function calib = loadCalibrationData(fn)

if nargin %if a filename is passed
    if isempty(fn)
        [fn, dr] = uigetfile('E:\SLAPmidata\Calibration\*.cal', 'Select a Calibration file');
        fn = [dr filesep fn];
    end
    calib = [];
    load(fn, '-mat')
else
    warning('Creating a DEFAULT calibration, for testing only! Use [] as input to load a user-selected file')
    sep = (8/sqrt(2))/2;  %10.8 or 8
    calib.galvos.offset.line1.X = sep;
    calib.galvos.offset.line1.Y = sep;
    calib.galvos.offset.line2.X = -sep;
    calib.galvos.offset.line2.Y = -sep;
    calib.galvos.offset.line3.X = sep;
    calib.galvos.offset.line3.Y = sep;
    calib.galvos.offset.line4.X = -sep;
    calib.galvos.offset.line4.Y = -sep;
    calib.galvos.offset.raster.X = 0;
    calib.galvos.offset.raster.Y = 0;
    calib.galvos.offset.E = [0 1];
    calib.galvos.offset.F = [0 1];
    calib.galvos.linelength = 4.3;  %100um is about 2.15 rad
    calib.galvos.VperDeg = 1;  %the voltage on the galvos that corresponds to 1 optical degree.
    calib.galvos.factor = [8 8 8 8 8 8]; %The number of V of galvo AI that corresponds to 1V of galvo AO; AI is smaller than AO.
    
    
    calib.pockels.fast = [0 1.6]; %Voltages for AB galvo / CD galvo respectively
    calib.pockels.slow = [0 1.6]; %Voltages for beam OFF/ON respectively
    
    calib.pockels.deadsamps = 0;

end