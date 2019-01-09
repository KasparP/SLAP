function SLAPMi_laserpower (scandata, hSLAPMi)

F = [scandata.frames.pmtData];

for line = 1:4
    I(line) = nanmean(nansum(F(scandata.line==line,:),1),2);
end

scaleby = min(I) ./ I;

disp('This data trace needs to be scaled by:')
disp(scaleby);
%change scaleby
scaleby = hSLAPMi.calib.B.scaleby .* scaleby;
scaleby = scaleby ./ min(scaleby);

disp('The new global scaling is:')
disp(scaleby)

if strcmpi( questdlg('Update calibration?'), 'Yes')
    hSLAPMi.calib.B.scaleby = scaleby;
    %ask whether to save calibration
    if strcmpi( questdlg('Save calibration?'), 'Yes')
        calib = hSLAPMi.calib;
        [fn, dr] = uiputfile([hSLAPMi.dataDir filesep 'calibration' filesep '*.cal']);
        save([dr filesep fn], 'calib');
    end
end
end