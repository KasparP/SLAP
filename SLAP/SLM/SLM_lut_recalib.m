function newlut=SLM_lut_recalib(lut,MLO_lut)
% recalibrate lut for MLO-provided lut such that ODP works
pix=MLO_lut(:,1);
volt=MLO_lut(:,2);
luton=lut(:,:,2);
lutoff=lut(:,:,1);

on=zeros(size(luton));
off=zeros(size(lutoff));
for i=1:length(on(:))
    [~,onInd]=min(abs(volt-luton(i)));
    [~,offInd]=min(abs(volt-lutoff(i)));
    on(i)=onInd-1;
    off(i)=offInd-1;
end
newlut=cat(3,off,on);
end