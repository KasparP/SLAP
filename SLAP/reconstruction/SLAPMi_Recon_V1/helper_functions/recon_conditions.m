function recon_conditions(sys_recon)
fprintf('Saving reconstruction conditions ... ');
appendname = AppendName(sys_recon);
% cd(sys_recon.opts.dr);
% fname = [sys_recon.opts.filename(1:end-11) 'RECON_' appendname];
% save(fname,'sys_recon','-v6')
mkdir(sys_recon.opts.dr, 'Recon');

fn = [sys_recon.opts.dr filesep 'Recon' filesep sys_recon.opts.filename(1:end-11) 'RECON_' appendname];
tag = [];
tagn = 0;
while exist([fn tag '.mat'], 'file')
    tagn = tagn+1;
    tag = ['_' int2str(tagn)];
end
save([fn tag '.mat'],'sys_recon','-v7.3')
fprintf('done \n');
end