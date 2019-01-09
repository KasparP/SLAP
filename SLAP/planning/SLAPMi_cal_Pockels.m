function calib = SLAPMi_cal_Pockels(calib)

%Some functions that populate the field calib.B, which deals with the
%brightness of illumination
curve = medfilt2(calib.B.curve(:), [5 1]);
curve = curve./max(curve);
[~,maxind] = max(curve);
V = 0:0.01:2;
calib.B.V2B = griddedInterpolant(V(1:maxind), curve(1:maxind));
calib.B.B2V = griddedInterpolant(curve(1:maxind), V(1:maxind));

end

