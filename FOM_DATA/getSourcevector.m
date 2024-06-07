%%
function f_1 = getSourcevector(HFM,var_para)
fvec = var_para.gx(HFM.P(:,1),HFM.P(:,2),var_para.x0s,var_para.y0s,var_para.sig0,var_para.sigx,var_para.sigy);
f_1 = zeros(HFM.DOF_full,1);
f_1(1:2:end,1) = fvec; f_1(2:2:end,1) = fvec;

f_1 = HFM.M_rho0*f_1(HFM.free_dofs,1);

end
