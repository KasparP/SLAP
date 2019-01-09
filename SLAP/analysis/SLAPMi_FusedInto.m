function FusedInfoDataset = SLAPMi_FusedInto
%select files
[fns, dr] = uigetfile('*PROBLEMDATA.mat', 'Select your saved RECON files for analysis', 'multiselect', 'on');
if ~iscell(fns)
    fns = {fns};
end
fns = sort_nat(fns);


for fnum = 1:length(fns)
    disp(['Loading Segmentation: ' int2str(fnum) ' of ' int2str(length(fns)) '...'])
    load([dr fns{fnum}],'S');
    if fnum==1
        FusedInfoDataset.fusedInto = nan(length(fnum),length(S.fusedInto));
    end
    FusedInfoDataset.fusedInto(fnum,:) = S.fusedInto;
end
mkdir([dr filesep 'SLAPMi_FusedInto']);
save([dr '/SLAPMi_FusedInto/FusedInto'],'FusedInfoDataset')

end