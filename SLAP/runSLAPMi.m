function runSLAPMi(sim)
    if nargin < 1
        sim = false;
    end
    
    cd('E:\SLAPmiData');
    
    if evalin('base','exist(''hSI'',''var'')')
        hSI = evalin('base','hSI');
        if ~most.idioms.isValidObj(hSI)
            hSI = [];
        end
    else
        hSI = [];
    end
    
    if evalin('base','exist(''hSLAPMi'',''var'')')
        hSLAPMi = evalin('base','hSLAPMi');
    else
        hSLAPMi = [];
    end
    
    if evalin('base','exist(''gSLAPMi'',''var'')')
        gSLAPMi = evalin('base','gSLAPMi');
    else
        gSLAPMi = [];
    end

    if most.idioms.isValidObj(hSLAPMi)
        if most.idioms.isValidObj(gSLAPMi)
            gSLAPMi.raise();
        else
            assignin('base', 'gSLAPMi', slapmiGui(hSLAPMi));
        end
    else
        hSLAPMi = slapmi(hSI,sim);
        assignin('base', 'hSLAPMi', hSLAPMi);
        assignin('base', 'gSLAPMi', slapmiGui(hSLAPMi));
    end
end
