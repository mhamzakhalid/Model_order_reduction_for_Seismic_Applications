function [wr,wi] =  getoptimalparameterWeeks(M_HFM,K_HFM,Vk,bi_HFM,smax,Rec,var_para)
sig0 = 0.12;
sigmax = 0.6;
var_para.si_max = smax;
bmax = var_para.si_max;
Ki = Vk'*(K_HFM*Vk);
Mi = Vk'*(M_HFM*Vk);
bi = Vk'*(bi_HFM);
Ntol = 40;

tic
[wr,wi] = wpar2([], var_para.t(end), var_para.Nz, sig0, sigmax, bmax, Ntol,Mi,Ki,bi,Vk,Rec,var_para);
toc

end