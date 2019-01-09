function sys = calc_corr( sys )
sys.recon_y = sys.output.PS*sys.output.F;
c = zeros(size(sys.input.y,1),1);
for i = 1:size(sys.input.y,1)
    cc = corrcoef(sys.input.y(i,:),sys.recon_y(i,:));
    c(i) = cc(2);
end
sys.control_params.corr_coeff = nanmean(c);

% sys.recon_y = normr(sys.PS*sys.output.F);
% sys.control_params.corr_coeff = nanmean(diag(sys.recon_y*normr(sys.input.y)'));

end

